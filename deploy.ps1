# Redeploy the site to Netlify.
# Copies all root HTML pages into site/ (index.html lives in site/ only), then publishes.
Copy-Item "$PSScriptRoot\*.html" "$PSScriptRoot\site\"
netlify deploy --prod --dir "$PSScriptRoot\site"
