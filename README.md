# seowony.github.io

Seo-Won Chang's personal website. A single static page — no build step, no dependencies, no
framework. Edit `index.html`, commit, push.

Live at **https://seowony.github.io/** · CV at **https://seowony.github.io/cv.pdf**

```
index.html              the whole site (HTML + CSS + inline SVG illustrations)
cv.pdf                  CV, compiled from CV_Seo-Won_Chang.tex
CV_Seo-Won_Chang.tex    CV source
assets/photo.jpg        portrait, 560 px wide
.nojekyll               tells GitHub Pages to serve the files as-is
```

## First push

From inside this folder:

```bash
git init -b main
git add .
git commit -m "Initial site"
git remote add origin https://github.com/seowony/seowony.github.io.git
git push -u origin main
```

If the repository does not exist yet, create it on GitHub first: **New repository**, owner `seowony`,
name exactly `seowony.github.io`, **Public**, and do not add a README or .gitignore (this folder
already has the files).

Then on GitHub: **Settings → Pages → Build and deployment → Source: Deploy from a branch**, branch
`main`, folder `/ (root)`, Save. The site is live a minute or two later.

If the repository already exists and has commits, replace the last two lines above with:

```bash
git remote add origin https://github.com/seowony/seowony.github.io.git
git pull --rebase origin main
git push -u origin main
```

## Later changes

```bash
git add -A && git commit -m "what changed" && git push
```

GitHub Pages rebuilds within a minute. A hard refresh (⌘⇧R) clears the browser cache.

## Updating the content

- **CV** — recompile and replace `cv.pdf`:

  ```bash
  pdflatex CV_Seo-Won_Chang.tex && mv CV_Seo-Won_Chang.pdf cv.pdf
  ```

- **Photo** — replace `assets/photo.jpg`. Keep it under ~200 KB; 500–700 px wide is plenty.
- **Text** — plain HTML blocks in `index.html`, each under a comment banner
  (`<!-- ===== RECENT ===== -->`, `<!-- ===== RUBIN / LSST ===== -->`, and so on). Copy an existing
  block and edit it.
- **Colours** — the values at the top of the `<style>` block (`--paper`, `--ink`, `--accent`,
  `--label`) control the whole page.

## Still to do

- "Pictures from the work" — paper figures (KS4 DR1 coverage, BLAP light curve, BOES flare spectra,
  the Rubin EDP2 → 7DT diagram) were left out of this version.
- Per-topic research pages, one for each of the four research cards.
- A custom domain, if wanted: add a `CNAME` file containing the domain, and point a DNS CNAME record
  at `seowony.github.io`.
