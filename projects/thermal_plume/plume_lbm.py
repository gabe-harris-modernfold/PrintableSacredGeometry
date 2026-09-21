"""Thermal lattice-Boltzmann solver for the natural-convection plume off a standing body.

D2Q9 for momentum with a two-relaxation-time (TRT) collision and Guo forcing; D2Q5 for
temperature. Buoyancy enters through the Boussinesq body force only -- there is no
radiative transfer in the fluid, because air is very nearly transparent at the 9-10 um
wavelengths skin emits. The measured radiative coefficients in body_model.py describe
energy leaving for the walls, and are accounted for separately in the writeup rather
than driving this flow.

Wall temperature is per-segment and comes from measurement (see body_model.SEGMENTS),
so the buoyancy source varies along the body and is deliberately non-monotonic: the
sole of the foot sits 5.6 K above ambient while the anterior neck sits 10.6 K above it.

Run it:
    python plume_lbm.py --ra 1e7 --seconds 60 --out run_ra1e7
    python plume_lbm.py --calibrate            # 8 s on a coarse grid, prints stability

Everything is nondimensional inside the solver. Lengths are lattice cells, velocity is
scaled on the buoyant free-fall velocity sqrt(g beta dT L), and temperature on the
neck's excess over ambient. Physical units are attached only on output.
"""

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import body_model as bm

# ---------------------------------------------------------------------------
# Lattice
# ---------------------------------------------------------------------------

CX = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1], dtype=np.int64)
CY = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1], dtype=np.int64)
W = np.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36])
OPP = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])

TCX = np.array([0, 1, 0, -1, 0], dtype=np.int64)
TCY = np.array([0, 0, 1, 0, -1], dtype=np.int64)
TW = np.array([1 / 3, 1 / 6, 1 / 6, 1 / 6, 1 / 6])
TOPP = np.array([0, 3, 4, 1, 2])

CS2 = 1.0 / 3.0


class PlumeSim:
    def __init__(self, ra=1.0e7, pr=0.71, u_lat=0.05, dx=0.010,
                 width_m=2.00, height_m=3.20, smagorinsky=0.10,
                 magic=3.0 / 16.0, seed=20260830, sponge_m=0.60):
        self.ra, self.pr, self.u_lat, self.dx = ra, pr, u_lat, dx
        self.smagorinsky, self.magic = smagorinsky, magic

        self.nx = int(round(width_m / dx))
        self.ny = int(round(height_m / dx))
        self.l_body = bm.BODY_HEIGHT / dx          # body height in cells

        # Transport coefficients from Ra and Pr, with the body height as the length.
        re = np.sqrt(ra / pr)
        self.nu = u_lat * self.l_body / re
        self.alpha = self.nu / pr
        self.tau_p = 3.0 * self.nu + 0.5
        self.tau_m = 0.5 + magic / (self.tau_p - 0.5)
        self.tau_g = 3.0 * self.alpha + 0.5
        self.g_lat = u_lat ** 2 / self.l_body      # since u_lat = sqrt(g beta dT L)

        # Physical scales.
        self.u_ff = bm.free_fall_velocity()
        self.dt = (u_lat / self.u_ff) * dx
        self.t_buoy = bm.BODY_HEIGHT / self.u_ff

        # Geometry: body midline centred, feet one cell above the floor row.
        self.body_x = self.nx // 2
        self.body_y = 1
        solid, seg_id, theta_wall = bm.rasterise(
            self.nx, self.ny, dx, self.body_x, self.body_y)
        solid[0, :] = True                          # the floor
        seg_id[0, :] = -1
        theta_wall[0, :] = 0.0                      # floor at ambient
        self.solid, self.seg_id, self.theta_wall = solid, seg_id, theta_wall
        self.fluid = ~solid

        # Only solid cells that actually touch fluid carry a wall condition. Applying
        # anti-bounce-back to the body's interior lets neighbouring solid cells reflect
        # off each other with eigenvalue -1, which grows an alternating-sign mode in
        # theta without bound. The interior is held at equilibrium instead.
        # A population leaving a solid node in direction d is reflected only if the
        # neighbour it would land on is fluid. Reflecting into another solid cell --
        # which happens in the gap between the feet and along the foot/floor
        # junction -- lets adjacent solid cells bounce off each other, and the
        # anti-bounce-back's eigenvalue of -1 then grows an alternating mode in
        # theta without bound. Everything else in the body is held at equilibrium.
        def nbr_is_fluid(cx, cy):
            return np.roll(np.roll(self.fluid, -int(cy), axis=0), -int(cx), axis=1)

        self.reflect_f = np.array([solid & nbr_is_fluid(CX[d], CY[d])
                                   for d in range(9)])
        self.reflect_g = np.array([solid & nbr_is_fluid(TCX[d], TCY[d])
                                   for d in range(5)])
        self.shell = self.reflect_f[1:].any(axis=0)

        # Flat index arrays for the hot loop. Boolean-mask assignment on a
        # (9, ny, nx) array is far slower than 1D fancy indexing on a ravelled view.
        self.rf_idx = [np.nonzero(self.reflect_f[d].ravel())[0] for d in range(9)]
        self.rg_idx = [np.nonzero(self.reflect_g[d].ravel())[0] for d in range(5)]
        self.solid_idx = np.nonzero(solid.ravel())[0]
        self.theta_wall_solid = theta_wall.ravel()[self.solid_idx]
        self.tw_g = [TW[d] * theta_wall.ravel()[self.rg_idx[d]] for d in range(5)]
        # TRT works on opposite pairs; the rest population is handled separately.
        self.pairs = [(1, 3), (2, 4), (5, 7), (6, 8)]

        self._build_wall_links()

        # Sponge layer below the open top. Without it the zero-gradient outflow
        # reflects, and the plume's peak velocity parks on the last row and keeps
        # climbing instead of settling. It sits well above every probe.
        self.sponge_cells = int(round(sponge_m / dx))
        sigma = np.zeros((self.ny, 1))
        if self.sponge_cells > 0:
            j0 = self.ny - self.sponge_cells
            ramp = (np.arange(j0, self.ny) - j0) / max(1, self.sponge_cells - 1)
            sigma[j0:, 0] = 0.35 * ramp ** 2
        self.sponge = sigma[None, :, :]
        self.sponge_top_m = (self.ny - self.sponge_cells - self.body_y) * dx

        rng = np.random.default_rng(seed)
        self.rho = np.ones((self.ny, self.nx))
        self.ux = np.zeros((self.ny, self.nx))
        self.uy = np.zeros((self.ny, self.nx))
        # A still room is a fiction; seed a small disturbance field instead.
        self.ux[self.fluid] = 1e-4 * rng.standard_normal(self.fluid.sum())
        self.uy[self.fluid] = 1e-4 * rng.standard_normal(self.fluid.sum())
        self.theta = np.zeros((self.ny, self.nx))

        self.f = self._feq(self.rho, self.ux, self.uy)
        self.g = self._geq(self.theta, self.ux, self.uy)
        self.step_count = 0

        # Probes.
        head_top = self.body_y + int(round(1.695 / dx))
        self.probe_row = min(head_top + int(round(0.50 / dx)), self.ny - 3)
        self.probe_col = self.body_x
        self.chest_row = self.body_y + int(round(1.25 / dx))
        self.ankle_row = self.body_y + int(round(0.12 / dx))
        self.history = {"t": [], "v_probe": [], "v_max": [], "ke": [],
                        "v_row_max": []}

    # -- lattice helpers ----------------------------------------------------

    def _feq(self, rho, ux, uy):
        cu = CX[:, None, None] * ux + CY[:, None, None] * uy
        usq = ux * ux + uy * uy
        return W[:, None, None] * rho * (
            1.0 + cu / CS2 + 0.5 * cu * cu / CS2 ** 2 - 0.5 * usq / CS2)

    def _geq(self, theta, ux, uy):
        cu = TCX[:, None, None] * ux + TCY[:, None, None] * uy
        return TW[:, None, None] * theta * (1.0 + cu / CS2)

    @staticmethod
    def _stream(h, cx, cy):
        out = np.empty_like(h)
        for i in range(h.shape[0]):
            out[i] = np.roll(np.roll(h[i], int(cy[i]), axis=0), int(cx[i]), axis=1)
        return out

    def _build_wall_links(self):
        """Fluid cells touching the body, tagged with the segment they touch.

        Used for the per-segment Nusselt estimate. Only the four axial links are
        counted, each contributing equally to the segment's exposed perimeter.
        """
        links = {i: [] for i in range(len(bm.SEGMENTS))}
        for d in range(1, 5):
            dy, dx_ = int(TCY[d]), int(TCX[d])
            nbr_seg = np.roll(np.roll(self.seg_id, -dy, axis=0), -dx_, axis=1)
            nbr_solid = np.roll(np.roll(self.solid, -dy, axis=0), -dx_, axis=1)
            exposed = self.fluid & nbr_solid & (nbr_seg >= 0)
            rows, cols = np.nonzero(exposed)
            for r, c in zip(rows, cols):
                links[int(nbr_seg[r, c])].append((r, c))
        self.wall_links = {k: (np.array([p[0] for p in v]), np.array([p[1] for p in v]))
                           for k, v in links.items() if v}

    # -- one step -----------------------------------------------------------

    def step(self):
        f, g = self.f, self.g

        rho = f[0].copy()
        for d in range(1, 9):
            rho += f[d]
        theta = g[0].copy()
        for d in range(1, 5):
            theta += g[d]
        fy = self.g_lat * theta * rho          # Boussinesq buoyancy, upward

        mx = f[1] - f[3] + f[5] - f[6] - f[7] + f[8]
        my = f[2] - f[4] + f[5] + f[6] - f[7] - f[8]
        ux = mx / rho
        uy = (my + 0.5 * fy) / rho
        ux[self.solid] = 0.0
        uy[self.solid] = 0.0
        usq = ux * ux + uy * uy

        tau_p = self.tau_p
        if self.smagorinsky > 0.0:
            feq = self._feq(rho, ux, uy)
            fneq = f - feq
            pxx = (fneq * (CX * CX)[:, None, None]).sum(axis=0)
            pyy = (fneq * (CY * CY)[:, None, None]).sum(axis=0)
            pxy = (fneq * (CX * CY)[:, None, None]).sum(axis=0)
            q = np.sqrt(pxx * pxx + pyy * pyy + 2.0 * pxy * pxy)
            tau_p = 0.5 * (self.tau_p + np.sqrt(
                self.tau_p ** 2 + 18.0 * np.sqrt(2.0) * self.smagorinsky ** 2 * q / rho))
        tau_m = 0.5 + self.magic / (tau_p - 0.5)
        inv_p, inv_m = 1.0 / tau_p, 1.0 / tau_m

        # Guo forcing. The prefactor takes tau_minus, not tau_plus: momentum is an odd
        # moment, so under TRT it relaxes at tau_minus, and only that choice makes the
        # per-step momentum increment come out equal to F. Using tau_plus here delivers
        # a few percent of the intended buoyancy and the plume never starts.
        gpre = 1.0 - 0.5 * inv_m

        f_post = np.empty_like(f)

        # Rest population: no antisymmetric part, and c_0 = 0 kills the cu terms.
        eq0 = W[0] * rho * (1.0 - 0.5 * usq / CS2)
        f_post[0] = f[0] - (f[0] - eq0) * inv_p + W[0] * gpre * (-uy / CS2) * fy

        # Opposite pairs. eq_plus and eq_minus have closed forms that avoid ever
        # materialising the full (9, ny, nx) equilibrium.
        for d, o in self.pairs:
            cu = CX[d] * ux + CY[d] * uy
            wr = W[d] * rho
            eq_plus = wr * (1.0 + 0.5 * cu * cu / CS2 ** 2 - 0.5 * usq / CS2)
            eq_minus = wr * cu / CS2
            f_plus = 0.5 * (f[d] + f[o])
            f_minus = 0.5 * (f[d] - f[o])
            relax_p = (f_plus - eq_plus) * inv_p
            relax_m = (f_minus - eq_minus) * inv_m
            common = W[d] * gpre * fy
            src_d = common * ((CY[d] - uy) / CS2 + cu * CY[d] / CS2 ** 2)
            src_o = common * ((-CY[d] - uy) / CS2 + cu * CY[d] / CS2 ** 2)
            f_post[d] = f[d] - relax_p - relax_m + src_d
            f_post[o] = f[o] - relax_p + relax_m + src_o

        g_post = np.empty_like(g)
        inv_g = 1.0 / self.tau_g
        for d in range(5):
            cu = TCX[d] * ux + TCY[d] * uy
            geq = TW[d] * theta * (1.0 + cu / CS2)
            g_post[d] = g[d] - (g[d] - geq) * inv_g

        # Hold the whole body at rest equilibrium, then overwrite only the populations
        # that stream out into fluid: halfway bounce-back for momentum, halfway
        # anti-bounce-back for temperature, which pins the wall at the segment's
        # measured skin temperature.
        si = self.solid_idx
        for d in range(9):
            fpd = f_post[d].reshape(-1)
            fpd[si] = W[d]
        for d in range(5):
            gpd = g_post[d].reshape(-1)
            gpd[si] = TW[d] * self.theta_wall_solid
        for d in range(1, 9):
            idx = self.rf_idx[d]
            if idx.size:
                f_post[d].reshape(-1)[idx] = f[OPP[d]].reshape(-1)[idx]
        for d in range(1, 5):
            idx = self.rg_idx[d]
            if idx.size:
                g_post[d].reshape(-1)[idx] = (-g[TOPP[d]].reshape(-1)[idx]
                                              + 2.0 * self.tw_g[d])

        f = self._stream(f_post, CX, CY)
        g = self._stream(g_post, TCX, TCY)

        if self.sponge_cells > 0:
            rho_s = f.sum(axis=0)
            f = (1.0 - self.sponge) * f + self.sponge * self._feq(
                rho_s, np.zeros_like(rho_s), np.zeros_like(rho_s))
            g = (1.0 - self.sponge) * g          # relaxes theta toward ambient

        # Open top and sides: zero gradient. Floor is solid, handled above.
        f[:, -1, :] = f[:, -2, :]
        g[:, -1, :] = g[:, -2, :]
        f[:, :, 0] = f[:, :, 1]
        f[:, :, -1] = f[:, :, -2]
        g[:, :, 0] = g[:, :, 1]
        g[:, :, -1] = g[:, :, -2]

        self.f, self.g = f, g
        self.rho, self.ux, self.uy, self.theta = rho, ux, uy, theta
        self.step_count += 1

        if self.step_count % 200 == 0:
            mean_rho = self.rho[self.fluid].mean()
            if np.isfinite(mean_rho) and mean_rho > 0:
                self.f *= 1.0 / mean_rho

    # -- diagnostics --------------------------------------------------------

    def stable(self):
        umag = np.hypot(self.ux, self.uy)
        return bool(np.isfinite(umag).all() and umag.max() < 0.4
                    and np.isfinite(self.theta).all() and np.abs(self.theta).max() < 5.0)

    def record(self):
        umag = np.hypot(self.ux, self.uy)
        self.history["t"].append(self.step_count * self.dt)
        self.history["v_probe"].append(float(self.uy[self.probe_row, self.probe_col]))
        self.history["v_max"].append(float(umag[self.fluid].max()))
        self.history["ke"].append(float((umag[self.fluid] ** 2).mean()))
        # Peak upward velocity across the probe row -- the like-for-like comparison
        # with the measured 0.2-0.3 m/s at 0.5 m above the head. The global maximum
        # sits in the near-wall jet and is a different quantity.
        self.history["v_row_max"].append(float(self.uy[self.probe_row].max()))

    def segment_nusselt(self):
        """Perimeter-averaged wall temperature gradient per segment.

        Reported as a dimensionless conductance, (theta_wall - theta_fluid)/theta_wall
        averaged over the segment's exposed links. Absolute values are not comparable
        to a 3D manikin's W/m^2K; the ratios between segments are the test.
        """
        out = {}
        for idx, (rows, cols) in self.wall_links.items():
            tw = (bm.SEGMENTS[idx]["T_skin"] - bm.T_AMBIENT) / bm.DELTA_T_REF
            tf = self.theta[rows, cols]
            out[bm.SEGMENTS[idx]["name"]] = float(np.mean((tw - tf) / tw))
        return out

    def boundary_layer_profile(self, threshold=0.10):
        """Half-thickness of the rising layer versus height, in metres.

        Measured on the body's right side as the distance from the silhouette edge
        out to where upward velocity falls below `threshold` of the local maximum.
        """
        heights, widths = [], []
        for j in range(self.body_y, min(self.ny - 1, self.body_y + int(1.70 / self.dx))):
            row_solid = np.nonzero(self.solid[j])[0]
            row_solid = row_solid[row_solid > self.body_x] if row_solid.size else row_solid
            if row_solid.size == 0:
                continue
            edge = int(row_solid.max())
            prof = self.uy[j, edge + 1:]
            if prof.size < 5 or prof.max() <= 0:
                continue
            cut = threshold * prof.max()
            beyond = np.nonzero(prof < cut)[0]
            if beyond.size == 0:
                continue
            heights.append((j - self.body_y) * self.dx)
            widths.append(float(beyond[0]) * self.dx)
        return np.array(heights), np.array(widths)

    def metadata(self):
        return dict(
            ra=self.ra, pr=self.pr, u_lat=self.u_lat, dx=self.dx, dt=self.dt,
            nx=self.nx, ny=self.ny, l_body_cells=self.l_body,
            nu_lat=self.nu, alpha_lat=self.alpha,
            tau_p=self.tau_p, tau_m=self.tau_m, tau_g=self.tau_g,
            g_lat=self.g_lat, smagorinsky=self.smagorinsky,
            sponge_cells=self.sponge_cells, sponge_top_m=self.sponge_top_m,
            u_ff=self.u_ff, t_buoy=self.t_buoy,
            ra_real=bm.rayleigh(bm.BODY_HEIGHT, bm.DELTA_T_REF),
            body_x=self.body_x, body_y=self.body_y,
            probe_row=self.probe_row, probe_col=self.probe_col,
            chest_row=self.chest_row, ankle_row=self.ankle_row,
            T_ambient=bm.T_AMBIENT, delta_T_ref=bm.DELTA_T_REF,
        )


def run(sim, seconds, out_dir, snapshot_every_s=0.25, quiet=False):
    os.makedirs(out_dir, exist_ok=True)
    total = int(round(seconds / sim.dt))
    snap_every = max(1, int(round(snapshot_every_s / sim.dt)))
    record_every = max(1, snap_every // 8)

    snaps_u, snaps_v, snaps_t, snaps_time = [], [], [], []
    nusselt_acc, nusselt_n = {}, 0
    t0 = time.time()
    diverged = None

    # A long run that blows up at minute 30 must still hand back the 30 minutes it
    # earned, so divergence ends the loop and is recorded rather than raised.
    for n in range(total):
        sim.step()
        if n % record_every == 0:
            sim.record()
            if not sim.stable() or not np.isfinite(sim.rho).all():
                diverged = dict(step=n, t=n * sim.dt,
                                umax=float(np.hypot(sim.ux, sim.uy).max()))
                print(f"  DIVERGED at step {n}, t={n * sim.dt:.2f}s "
                      f"({100 * n / total:.0f}% of the run) -- saving what exists",
                      flush=True)
                break
        if n % snap_every == 0:
            snaps_u.append(sim.ux.astype(np.float32))
            snaps_v.append(sim.uy.astype(np.float32))
            snaps_t.append(sim.theta.astype(np.float32))
            snaps_time.append(n * sim.dt)
            # Average the wall conductance over the second half only.
            if n * sim.dt > 0.5 * seconds:
                for k, v in sim.segment_nusselt().items():
                    nusselt_acc[k] = nusselt_acc.get(k, 0.0) + v
                nusselt_n += 1
        if not quiet and n % max(1, total // 20) == 0 and n:
            el = time.time() - t0
            print(f"  {100 * n / total:5.1f}%  t={n * sim.dt:6.2f}s  "
                  f"|u|max={np.hypot(sim.ux, sim.uy).max():.4f}  "
                  f"{n / el:6.0f} steps/s  eta {el * (total - n) / n / 60:.1f} min",
                  flush=True)

    bl_h, bl_w = sim.boundary_layer_profile()
    nusselt = {k: v / nusselt_n for k, v in nusselt_acc.items()} if nusselt_n else {}

    np.savez_compressed(
        os.path.join(out_dir, "fields.npz"),
        u=np.array(snaps_u), v=np.array(snaps_v), theta=np.array(snaps_t),
        snap_t=np.array(snaps_time), solid=sim.solid, seg_id=sim.seg_id,
        hist_t=np.array(sim.history["t"]),
        hist_vprobe=np.array(sim.history["v_probe"]),
        hist_vmax=np.array(sim.history["v_max"]),
        hist_ke=np.array(sim.history["ke"]),
        hist_vrowmax=np.array(sim.history["v_row_max"]),
        bl_height=bl_h, bl_width=bl_w)

    meta = sim.metadata()
    meta.update(seconds=seconds, steps=total, wall_clock_s=time.time() - t0,
                segment_conductance=nusselt, snapshots=len(snaps_time),
                diverged=diverged,
                seconds_completed=float(snaps_time[-1]) if snaps_time else 0.0)
    with open(os.path.join(out_dir, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    return meta


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ra", type=float, default=1.0e7)
    ap.add_argument("--seconds", type=float, default=60.0)
    ap.add_argument("--dx", type=float, default=0.010)
    ap.add_argument("--u", type=float, default=0.05)
    ap.add_argument("--smag", type=float, default=0.10)
    ap.add_argument("--out", default="run")
    ap.add_argument("--calibrate", action="store_true",
                    help="short coarse run that only reports stability")
    args = ap.parse_args()

    if args.calibrate:
        for ra in (1.0e6, 3.0e6, 1.0e7, 3.0e7):
            sim = PlumeSim(ra=ra, dx=0.020, u_lat=args.u, smagorinsky=args.smag)
            print(f"Ra={ra:.0e}  grid {sim.nx}x{sim.ny}  L_body={sim.l_body:.0f} cells  "
                  f"nu={sim.nu:.2e}  tau_p={sim.tau_p:.4f}  tau_g={sim.tau_g:.4f}  "
                  f"dt={sim.dt * 1e3:.3f} ms")
            try:
                steps = int(round(8.0 / sim.dt))
                t0 = time.time()
                for n in range(steps):
                    sim.step()
                    if n % 400 == 0 and not sim.stable():
                        raise RuntimeError(f"unstable at t={n * sim.dt:.2f}s")
                rate = steps / (time.time() - t0)
                umax = np.hypot(sim.ux, sim.uy).max()
                print(f"          stable to 8 s. |u|max={umax:.4f} "
                      f"({umax / args.u:.2f} U_ff), {rate:.0f} steps/s")
            except RuntimeError as exc:
                print(f"          FAILED: {exc}")
        return

    sim = PlumeSim(ra=args.ra, dx=args.dx, u_lat=args.u, smagorinsky=args.smag)
    print(f"grid {sim.nx}x{sim.ny}, body {sim.l_body:.0f} cells, "
          f"tau_p={sim.tau_p:.4f}, tau_g={sim.tau_g:.4f}, dt={sim.dt * 1e3:.3f} ms")
    print(f"Ra_sim={args.ra:.1e} vs Ra_real={bm.rayleigh(bm.BODY_HEIGHT, bm.DELTA_T_REF):.2e} "
          f"(ratio {bm.rayleigh(bm.BODY_HEIGHT, bm.DELTA_T_REF) / args.ra:.0f})")
    meta = run(sim, args.seconds, args.out)
    print(f"done in {meta['wall_clock_s'] / 60:.1f} min -> {args.out}/")


if __name__ == "__main__":
    main()
