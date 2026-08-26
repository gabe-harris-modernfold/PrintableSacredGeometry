import numpy as np, trimesh, pyvista as pv, math
import spiral as S, params as P

def pvm(t):
    f = np.hstack([np.full((len(t.faces),1),3), t.faces]).ravel()
    return pv.PolyData(t.vertices.copy(), f)

sp = trimesh.load("drip_gutters.stl")
tap = trimesh.load("drip_splitter.stl")

# ghost vessel: ovoid crown over a tapered 12-sided body
z = np.linspace(0, S.VESSEL_H, 200); r = S.vessel_r(z)
ang = np.linspace(0, 2*math.pi, 12, endpoint=False)
hexf = math.cos(math.pi/12)/np.cos((ang % (math.pi/6)) - math.pi/12)
V, F = [], []
for i, (zz, rr) in enumerate(zip(z, r)):
    for a, h in zip(ang, hexf):
        V.append([rr*h*math.cos(a), rr*h*math.sin(a), zz])
n = len(ang)
for i in range(len(z)-1):
    for j in range(n):
        a=i*n+j; b=i*n+(j+1)%n; c=(i+1)*n+(j+1)%n; d=(i+1)*n+j
        F += [[a,b,c],[a,c,d]]
ves = trimesh.Trimesh(np.array(V), np.array(F), process=False)

pl = pv.Plotter(off_screen=True, window_size=(1750,1100), shape=(1,3))
pl.subplot(0,0)
pl.add_mesh(pvm(ves), color="#dfeef5", opacity=.28, smooth_shading=True)
pl.add_mesh(pvm(sp), color="#3f9ad0", smooth_shading=True, specular=.7)
pl.add_mesh(pvm(tap), color="#c0662f", smooth_shading=True, specular=.7)
pl.add_text("vessel + 16 drip spirals + condenser tap", font_size=11)
pl.camera_position=[(1500,-1350,1250),(0,0,275),(0,0,1)]; pl.camera.zoom(1.15)

pl.subplot(0,1)
# clipped, with the crown ghosted: the point of this part is that the cone's seat
# TOUCHES the dome's inner face, and that reads nowhere else
crown = pvm(trimesh.load("vessel_shell_2.stl")).clip(normal="y", origin=(0,0,0))
pl.add_mesh(crown, color="#dfeef5", opacity=.40, smooth_shading=True)
pl.add_mesh(pvm(tap).clip(normal="y", origin=(0,0,0)),
            color="#c0662f", smooth_shading=True, specular=.7)
_zs = S.crown_in_z(S.SEAT_R); _th = np.linspace(math.pi, 2*math.pi, 80)
pl.add_mesh(pv.MultipleLines(np.c_[S.SEAT_R*np.cos(_th), S.SEAT_R*np.sin(_th),
                                   np.full(80, _zs)]), color="#18c46a", line_width=6)
pl.add_text("condenser tap: 10mm intake -> 2 -> 4 -> 8 -> 16 x 2.5mm" + chr(10)
            + f"seat touches the dome on r={S.SEAT_R:.0f} at z={_zs:.2f}",
            font_size=10)
pl.camera_position=[(300,-260,610),(0,0,524),(0,0,1)]; pl.camera.zoom(1.15)

pl.subplot(0,2)
g,_,_ = S.gutter(S.helix(0)[0])
pl.add_mesh(pvm(g), color="#3f9ad0", smooth_shading=True, specular=.6)
pl.add_text("one gutter: level runs, 10mm drops, notched lip", font_size=11)
pl.camera_position=[(560,-330,470),(60,-40,300),(0,0,1)]; pl.camera.zoom(1.6)
pl.screenshot("spiral_views.png")
print("ok")
