#!/usr/bin/env python3
"""
Build realbogus.html: a small "Real or bogus?" quiz from labelled IMSNG/GECKO candidates.

Input: a folder (default realbogus_input/, not tracked by git) with labels.csv and the stamp images.

    id,label,new,ref,diff,score,note,facility
    c001,real,c001_new.fits,c001_ref.fits,c001_diff.fits,0.97,SN 2025abc in NGC 1234,KCT
    c002,bogus,c002_new.png,c002_ref.png,c002_diff.png,0.08,Cosmic ray,RASA36

  label     real or bogus (the team's vetted label)
  new/ref/diff  image files: FITS (scaled with zscale) or PNG/JPG
  score     optional machine-learning real/bogus score (0-1)
  note      optional short explanation shown after answering
  facility  optional telescope name

Output: assets/realbogus/<id>.jpg (new | ref | diff side by side) and realbogus.html from
realbogus_template.html. Only use candidates that are public (e.g. reported to TNS) or plainly bogus.

    pip install astropy numpy pillow
    python tools/make_realbogus.py [--input realbogus_input] [--size 160]
"""
import argparse, csv, html, json, os, sys

import numpy as np
from PIL import Image

OUT = "assets/realbogus"


def load_stamp(path, size):
    if path.lower().endswith((".fits", ".fit", ".fz", ".fits.gz")):
        from astropy.io import fits
        from astropy.visualization import ZScaleInterval
        with fits.open(path) as h:
            data = next(x.data for x in h if getattr(x, "data", None) is not None and x.data.ndim == 2)
        data = np.nan_to_num(np.asarray(data, dtype=float))
        lo, hi = ZScaleInterval().get_limits(data)
        img = np.clip((data - lo) / max(hi - lo, 1e-9), 0, 1)
        im = Image.fromarray((img[::-1] * 255).astype(np.uint8))      # FITS row 0 is the bottom
    else:
        im = Image.open(path).convert("L")
    return im.resize((size, size), Image.LANCZOS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="realbogus_input")
    ap.add_argument("--size", type=int, default=160, help="pixels per stamp")
    args = ap.parse_args()
    if not os.path.exists("realbogus_template.html"):
        sys.exit("Run this from the website root (realbogus_template.html not found).")
    table = os.path.join(args.input, "labels.csv")
    if not os.path.exists(table):
        sys.exit(f"{table} not found")
    os.makedirs(OUT, exist_ok=True)

    items = []
    with open(table) as f:
        for row in csv.DictReader(f):
            cid, label = row["id"].strip(), row["label"].strip().lower()
            if label not in ("real", "bogus"):
                print(f"{cid}: label must be real or bogus, skipped")
                continue
            try:
                stamps = [load_stamp(os.path.join(args.input, row[k].strip()), args.size) for k in ("new", "ref", "diff")]
            except Exception as e:
                print(f"{cid}: {e}, skipped")
                continue
            gap = max(2, args.size // 40)
            strip = Image.new("L", (3 * args.size + 2 * gap, args.size), 255)
            for i, s in enumerate(stamps):
                strip.paste(s, (i * (args.size + gap), 0))
            strip.save(os.path.join(OUT, f"{cid}.jpg"), quality=88)
            score = row.get("score", "").strip()
            items.append(dict(id=cid, label=label, img=f"{OUT}/{cid}.jpg",
                              score=float(score) if score else None,
                              note=html.escape(row.get("note", "").strip()),
                              facility=html.escape(row.get("facility", "").strip())))
    if not items:
        sys.exit("No candidates: realbogus.html was not written.")
    n_real = sum(i["label"] == "real" for i in items)
    page = open("realbogus_template.html").read()
    page = page.replace("/* RB_DATA */", "const RB = " + json.dumps(items, separators=(",", ":")) + ";")
    page = page.replace("<!-- RB_COUNT -->", f"{len(items)} candidates ({n_real} real, {len(items) - n_real} bogus)")
    open("realbogus.html", "w").write(page)
    print(f"Wrote realbogus.html with {len(items)} candidates ({n_real} real).")


if __name__ == "__main__":
    main()
