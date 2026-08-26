

## Verification log

Three defects found by checking built geometry against spec, and fixed:

1. **drip_gutters winding.** Watertight but winding-inconsistent - 104 bad edge pairs
   (= 2 x 52 section points) at the cap/side junction. The C-section came out clockwise
   while the sweep's side quads and end-cap fans both assume counter-clockwise. Section
   is now forced CCW, with fix_winding/fix_normals as a safety net.
2. **Gutters through the shell.** Placed off the vessel's circumradius, but a 12-sided
   shell's inner face sits at circumradius x cos(pi/12) = 0.966x. Off by 3.4 mm on the
   body. Now placed off the inscribed radius; worst clearance +5.7 mm, per z-slice.
3. **Siphon height.** Standpipe originally topped out above the bed rim, so the tray
   would overflow before the siphon tripped. Now h=50: flood level 36 mm below the rim,
   bell cap 6 mm clear.

Two further 'failures' were artifacts of the checks themselves, not the design -
comparing max/min radii across a whole band rather than per z-slice, and measuring the
bed tray's overall bounds (which include its rim) instead of the siphon alone.
