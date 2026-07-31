# Regenerate menu-page preview images in site/previews/.
# Screenshots each root HTML page with headless Edge, then resizes to 800px WebP.
# Usage: .\make-previews.ps1            (all pages)
#        .\make-previews.ps1 -Only bifilar-rodin-coil-simulator   (one page)
param([string]$Only = "")

$edge = @("C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
          "C:\Program Files\Microsoft\Edge\Application\msedge.exe") | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $edge) { throw "Microsoft Edge not found." }

$tmp = Join-Path $env:TEMP "psg-previews"
New-Item -ItemType Directory -Force $tmp | Out-Null
New-Item -ItemType Directory -Force "$PSScriptRoot\site\previews" | Out-Null

$pages = Get-ChildItem "$PSScriptRoot\*.html"
if ($Only) { $pages = $pages | Where-Object { $_.BaseName -eq $Only } }

foreach ($page in $pages) {
    $name = $page.BaseName
    $png = "$tmp\$name.png"
    $fileUrl = "file:///" + $page.FullName.Replace('\', '/')
    & $edge --headless=new --disable-gpu --hide-scrollbars --window-size=1280,800 `
        --virtual-time-budget=12000 --screenshot="$png" $fileUrl 2>$null | Out-Null
    python -c @"
from PIL import Image
img = Image.open(r'$png').convert('RGB')
img = img.resize((800, round(img.height * 800 / img.width)), Image.LANCZOS)
out = r'$PSScriptRoot\site\previews\$name.webp'
img.save(out, 'WEBP', quality=82, method=6)
import os; print(f'{out}  {os.path.getsize(out)//1024} KB')
"@
}
Write-Output "Done. Remember to add new pages to site/index.html, then run deploy.ps1."
