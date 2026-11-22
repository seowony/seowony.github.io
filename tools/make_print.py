#!/usr/bin/env python3
"""
Make print files from a KS4 atlas: images in list order with their names.

  <atlas>_A4.pdf       A4 pages, 4 x 5 images per page (300 dpi)
  <atlas>_poster.pdf   one poster page (A0 at 300 dpi by default; prints well scaled to A1, A2 or A3)
  <atlas>_poster.jpg   the same poster as a JPEG

    python tools/make_print.py ngc --out print/
    python tools/make_print.py peculiar --poster A1 --out print/
    python tools/make_print.py ngc --exclude ngc2822,ngc5128 --out print/

Run it from the website root after tools/make_atlas.py (it reads assets/<atlas>/ and tools/lists/<atlas>.csv).
"""
import argparse, csv, math, os, sys

from PIL import Image, ImageDraw, ImageFont

DPI = 300
PAPER_MM = {"A0": (841, 1189), "A1": (594, 841), "A2": (420, 594), "A3": (297, 420), "A4": (210, 297)}
TITLES = {"ngc": "NGC galaxies in KS4 DR1", "peculiar": "Southern peculiar galaxies in KS4 DR1",
          "clusters": "Nearby galaxy clusters in KS4 DR1"}
CREDIT = ("KMTNet Synoptic Survey of Southern Sky (KS4) Data Release 1, BVI colour (Chang et al. 2026, JKAS). "
          "North up, east left. White bar: 1 arcmin.")
FONT_DIRS = ["/usr/share/fonts/truetype/dejavu/", "/Library/Fonts/", "/System/Library/Fonts/Supplemental/",
             "C:/Windows/Fonts/"]


def font(size, bold=False):
    names = (["DejaVuSans-Bold.ttf", "Arial Bold.ttf", "arialbd.ttf"] if bold else
             ["DejaVuSans.ttf", "Arial.ttf", "arial.ttf"])
    for d in FONT_DIRS:
        for n in names:
            if os.path.exists(d + n):
                return ImageFont.truetype(d + n, size)
    return ImageFont.load_default(size)


def px(mm):
    return int(round(mm / 25.4 * DPI))


def load_items(atlas):
    out_dir = f"assets/{atlas}"
    items = []
    with open(f"tools/lists/{atlas}.csv") as f:
        for row in csv.DictReader(f):
            path = os.path.join(out_dir, f"{row['id']}.jpg")
            if os.path.exists(path):
                items.append((row["name"], path))
    return items


def draw_grid(page, items, box, cols, rows, gap, cap_px, fnt):
    x0, y0, w, h = box
    cell = min((w - (cols - 1) * gap) / cols, (h - (rows - 1) * gap) / rows - cap_px)
    cell = int(cell)
    grid_w = cols * cell + (cols - 1) * gap
    xoff = x0 + (w - grid_w) // 2
    d = ImageDraw.Draw(page)
    for i, (name, path) in enumerate(items):
        r, c = divmod(i, cols)
        x = xoff + c * (cell + gap)
        y = y0 + r * (cell + cap_px + gap)
        im = Image.open(path).convert("RGB").resize((cell, cell), Image.LANCZOS)
        page.paste(im, (x, y))
        tw = d.textlength(name, font=fnt)
        d.text((x + (cell - tw) / 2, y + cell + cap_px * 0.18), name, fill=(30, 30, 30), font=fnt)


def a4_pages(atlas, items, cols=4, rows=5):
    W, H = px(210), px(297)
    m = px(12)
    per = cols * rows
    pages = []
    n_pages = math.ceil(len(items) / per)
    f_title, f_cap, f_foot = font(46, True), font(34), font(26)
    for p in range(n_pages):
        page = Image.new("RGB", (W, H), "white")
        d = ImageDraw.Draw(page)
        d.text((m, m), TITLES.get(atlas, atlas), fill=(20, 20, 20), font=f_title)
        chunk = items[p * per:(p + 1) * per]
        rng = f"{chunk[0][0]} to {chunk[-1][0]}"
        d.text((W - m - d.textlength(rng, font=f_cap), m + 8), rng, fill=(110, 110, 110), font=f_cap)
        top = m + px(10)
        bottom = H - m - px(9)
        draw_grid(page, chunk, (m, top, W - 2 * m, bottom - top), cols, rows, px(3), px(6), f_cap)
        d.text((m, H - m - 30), CREDIT, fill=(110, 110, 110), font=f_foot)
        pg = f"{p + 1} / {n_pages}"
        d.text((W - m - d.textlength(pg, font=f_foot), H - m - 30), pg, fill=(110, 110, 110), font=f_foot)
        pages.append(page)
    return pages


def poster(atlas, items, paper):
    wmm, hmm = PAPER_MM[paper]
    W, H = px(wmm), px(hmm)
    s = W / px(841)                                  # scale text with the paper size
    m = px(wmm * 0.03)
    f_title, f_sub, f_cap = font(int(150 * s), True), font(int(60 * s)), font(int(44 * s))
    head = int(330 * s)
    foot = int(110 * s)
    gap = int(30 * s)
    cap = int(70 * s)
    box = (m, m + head, W - 2 * m, H - 2 * m - head - foot)
    # pick the column count that gives the largest images
    best = None
    for cols in range(4, 30):
        rows = math.ceil(len(items) / cols)
        cell = min((box[2] - (cols - 1) * gap) / cols, (box[3] - (rows - 1) * gap) / rows - cap)
        if best is None or cell > best[0]:
            best = (cell, cols, rows)
    _, cols, rows = best
    page = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(page)
    d.text((m, m), TITLES.get(atlas, atlas), fill=(20, 20, 20), font=f_title)
    d.text((m, m + int(185 * s)), f"{len(items)} objects, {items[0][0]} to {items[-1][0]}", fill=(90, 90, 90), font=f_sub)
    draw_grid(page, items, box, cols, rows, gap, cap, f_cap)
    d.text((m, H - m - int(60 * s)), CREDIT, fill=(90, 90, 90), font=f_sub)
    return page, cols, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("atlas", choices=list(TITLES))
    ap.add_argument("--out", default="print")
    ap.add_argument("--poster", default="A0", choices=list(PAPER_MM))
    ap.add_argument("--exclude", default="", help="ids to leave out, e.g. ngc2822,ngc5128")
    args = ap.parse_args()
    skip = {x.strip().lower() for x in args.exclude.split(",") if x.strip()}
    items = [it for it in load_items(args.atlas)
             if os.path.basename(it[1])[:-4] not in skip]
    if not items:
        sys.exit(f"No images in assets/{args.atlas}/")
    os.makedirs(args.out, exist_ok=True)

    pages = a4_pages(args.atlas, items)
    a4 = os.path.join(args.out, f"{args.atlas}_A4.pdf")
    pages[0].save(a4, save_all=True, append_images=pages[1:], resolution=DPI, quality=92)
    print(f"{a4}: {len(pages)} pages, {len(items)} images")

    pst, cols, rows = poster(args.atlas, items, args.poster)
    pdf = os.path.join(args.out, f"{args.atlas}_poster.pdf")
    pst.save(pdf, resolution=DPI, quality=92)
    pst.save(os.path.join(args.out, f"{args.atlas}_poster.jpg"), quality=92, dpi=(DPI, DPI))
    print(f"{pdf}: {args.poster}, {cols} x {rows} grid, {pst.size[0]} x {pst.size[1]} px")


if __name__ == "__main__":
    main()
