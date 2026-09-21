# seowony.github.io

Personal website for Seo-Won Chang — https://seowony.github.io/

One static page. No build step, no dependencies. Edit `index.html`, commit, push.

```
index.html        the whole site (HTML, CSS and inline SVG in one file)
cv.tex / cv.pdf   web version of the CV — no referee details, no address, no grant amounts
assets/photo.jpg  portrait
```

Updating:

```bash
git add -A && git commit -m "what changed" && git push
```

To rebuild the CV after editing `cv.tex`: `pdflatex cv.tex`.

The full CV — with references, office address and grant amounts — is kept offline and is not
in this repository.
