"""Preview renders of the page, taken in headless Edge.

Each scene appends a few lines of script to a copy of web/language-of-the-plume.html that set up
the smoke, then screenshots it into previews/. Edge gets a throwaway profile directory; without
one, headless Edge hangs when a normal Edge window is already open.

    cd projects/language_of_the_plume && python shots.py            # every scene
    python shots.py air_day fork                                    # just these
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PAGE = os.path.join(HERE, "..", "..", "web", "language-of-the-plume.html")
EDGE = r"C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"
TAIL = "requestAnimationFrame(loop);\n</script>"

def adv(sec):
    return f"for(let i=0;i<60*{sec};i++){{sim.step(1/60);stepAir(sim,1/60);flowTick();}}"

def air(env, seed, sec=40):
    # hold off the communicator so the table shows the ridge's own air
    return (f"app.flow.env='{env}';sim.reset({{env:'{env}',gm:1,p:0.4,seed:{seed},k:WORLD.k}});"
            f"app.flow.next=1e9;" + adv(sec))

def stage(sec, signs="[]", env="dusk", seed=11, gm=0.6):
    return (f"setMode('free');sim.reset({{env:'{env}',gm:{gm},p:0.35,seed:{seed},k:0.65}});"
            f"sim.schedule({signs});for(let i=0;i<60*{sec};i++)sim.step(1/60);sim.measure();"
            "sim.step=()=>{};")  # freeze the frame so the screenshot catches this moment

# scene: (setup script, theme)
SCENES = {
    "sitting": ("", "dark"),
    "picked": ("app.flow.next=sim.t;" + adv(7) + "pick(app.flow.ep?app.flow.ep.truth:app.flow.opts[0],null);setConf(5);", "dark"),
    "air_day": (air("day", 9) + "pick(F.STONE,{dir:1,ring:3});", "dark"),
    "air_night": (air("night", 4), "dark"),
    "column": (stage(14, seed=3), "dark"),
    "tower": (stage(5.6, "[{form:F.TOWER,t0:0.5,T:7}]"), "dark"),
    "fork": (stage(6.5, "[{form:F.FORK,t0:3,T:6}]", seed=5, gm=0.8), "dark"),
    "ring": (stage(4.6, "[{form:F.RING,t0:1,T:6}]"), "dark"),
    "light": (stage(7.5, "[{form:F.TOWER,t0:0.5,T:7}]"), "light"),
}

def main(names):
    html = open(PAGE, encoding="utf-8").read()
    os.makedirs("previews", exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        for name in names:
            setup, theme = SCENES[name]
            page = html.replace(TAIL, "requestAnimationFrame(loop);\n" + setup + "\n</script>", 1)
            page = page.replace('<html lang="en">', f'<html lang="en" data-theme="{theme}">', 1)
            src = os.path.join(tmp, name + ".html")
            open(src, "w", encoding="utf-8").write(page)
            out = os.path.abspath(os.path.join("previews", name + ".png"))
            subprocess.run([EDGE, "--headless", "--user-data-dir=" + os.path.join(tmp, "profile"), "--disable-gpu",
                            "--hide-scrollbars", "--window-size=1200,900", "--virtual-time-budget=1500",
                            "--screenshot=" + out, "file:///" + src.replace("\\", "/")],
                           capture_output=True, timeout=90)
            print(name, "->", os.path.relpath(out) if os.path.exists(out) else "FAILED")

if __name__ == "__main__":
    main(sys.argv[1:] or list(SCENES))
