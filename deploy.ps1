# Redeploy the site to Netlify.
# Copies all web/ HTML pages into site/ (index.html lives in site/ only), then publishes.
Copy-Item "$PSScriptRoot\web\*.html" "$PSScriptRoot\site\"
netlify deploy --prod --dir "$PSScriptRoot\site"
