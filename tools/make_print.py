#!/usr/bin/env python3
"""
Make print files from a KS4 atlas: images in list order with their names.

  <atlas>_A4.pdf       A4 pages, 4 x 5 images per page (300 dpi)
  <atlas>_poster.pdf   one poster page (A0 at 300 dpi by default; prints well scaled to A1, A2 or A3)
  <atlas>_poster.jpg   the same poster as a JPEG (with --poster-jpg)
  <atlas>_tiles_3x2.pdf  (with --tiles 3x2) one poster split over A4 sheets that are pinned edge to edge:
                       every sheet holds whole images and the margins equal half a gap, so the grid runs on
                       across the joins. Print at 100% (actual size), no "fit to page".

    python tools/make_print.py ngc --out print/
    python tools/make_print.py peculiar --poster A1 --out print/
    python tools/make_print.py ngc --exclude ngc2822,ngc5128 --out print/

Run it from the website root after tools/make_atlas.py (it reads assets/<atlas>/ and tools/lists/<atlas>.csv).
"""
import argparse, csv, math, os, sys

from PIL import Image, ImageDraw, ImageFont

DPI = 300
PAPER_MM = {"A0": (841, 1189), "A1": (594, 841), "A2": (420, 594), "A3": (297, 420), "A4": (210, 297)}
TITLES = {"ngc": "NGC galaxies in KS4 DR1", "messier-caldwell": "Messier and Caldwell objects in KS4 DR1",
          "planetary": "Planetary nebulae in KS4 DR1",
          "nebulae": "Nebulae in KS4 DR1", "peculiar": "Southern peculiar galaxies in KS4 DR1",
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


def caption(d, name, x, y, width, fnt):
    """Centred caption under an image; long names ("C80 (NGC 5139, Omega Centauri)") go on two smaller lines."""
    if d.textlength(name, font=fnt) <= width:
        lines, f = [name], fnt
    else:
        if " (" in name:
            a, b = name.split(" (", 1)
            b = b.rstrip(")")
            lines = [a, b]
            if d.textlength(b, font=fnt.font_variant(size=int(fnt.size * 0.85))) > width and ", " in b:
                lines = [a] + b.split(", ", 1)        # number / designation / common name
        else:
            lines = [name]
        size = int(fnt.size * 0.85)
        f = fnt.font_variant(size=size)
        while size > fnt.size * 0.6 and max(d.textlength(t, font=f) for t in lines) > width:
            size -= 2
            f = fnt.font_variant(size=size)
    for i, t in enumerate(lines):
        tw = d.textlength(t, font=f)
        d.text((x + (width - tw) / 2, y + i * f.size * 1.1), t, fill=(30, 30, 30), font=f)


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
        caption(d, name, x, y + cell + cap_px * 0.18, cell, fnt)


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


def tiles(atlas, items, nx, ny, cell_mm=None, cap_mm=3.0):
    """Poster over nx x ny A4 portrait sheets with a seamless grid (each sheet has whole images)."""
    SW, SH = 210.0, 297.0
    # images per sheet: as many as fit with a cell of about 29 mm and margins of at least 5 mm
    cols, rows = 5, 7
    if cell_mm is None:
        cell_mm = min(SW / cols - 2 * 5.0, SH / rows - cap_mm - 2 * 5.0)
    gh = SW / cols - cell_mm                 # horizontal pitch minus cell; sheet margin = gh / 2
    gv = SH / rows - cell_mm - cap_mm
    W, H = px(SW), px(SH)
    total_cols, total_rows = cols * nx, rows * ny
    slots = total_cols * total_rows
    spare = slots - len(items)
    if spare < cols:
        sys.exit(f"{len(items)} images leave no room for the title in {slots} slots; use more sheets or --exclude")
    t_rows = max(1, min(3, spare // cols))           # title block: top rows of the top-left sheet
    n_credit = min(spare - t_rows * cols, cols * rows)  # credit block: last slots of the bottom-right sheet
    sheets = [[Image.new("RGB", (W, H), "white") for _ in range(nx)] for _ in range(ny)]
    f_cap, f_title, f_sub, f_tiny = font(32), font(92, True), font(40), font(20)
    cpx = px(cell_mm)

    def slot_xy(R, C):
        sy, sx = R // rows, C // cols
        lr, lc = R % rows, C % cols
        x = px(gh / 2 + lc * (cell_mm + gh))
        y = px(gv / 2 + lr * (cell_mm + cap_mm + gv))
        return sheets[sy][sx], x, y

    reserved = {(R, C) for R in range(t_rows) for C in range(cols)}
    last_sheet = [(R, C) for R in range((ny - 1) * rows, ny * rows) for C in range((nx - 1) * cols, nx * cols)]
    credit = last_sheet[len(last_sheet) - n_credit:] if n_credit else []
    reserved |= set(credit)
    free = [(R, C) for R in range(total_rows) for C in range(total_cols) if (R, C) not in reserved]

    page, x, y = slot_xy(0, 0)
    d = ImageDraw.Draw(page)
    block_w = W - 2 * px(gh / 2)
    size = 150
    while size > 60 and d.textlength(TITLES.get(atlas, atlas), font=font(size, True)) > block_w:
        size -= 6
    f_title, f_sub = font(size, True), font(int(size * 0.42))
    lines = [(TITLES.get(atlas, atlas), f_title, (20, 20, 20)),
             (f"{len(items)} {'galaxies' if atlas == 'ngc' else 'objects'}, {items[0][0]} to {items[-1][0]}", f_sub, (80, 80, 80)),
             ("KMTNet Synoptic Survey of Southern Sky (KS4) Data Release 1", f_sub, (80, 80, 80)),
             ("BVI colour. North up, east left. White bar: 1 arcmin.", f_sub, (80, 80, 80))]
    yy = y + px(6)
    for txt, f, col in lines:
        d.text((x, yy), txt, fill=col, font=f)
        yy += int(size * 1.45) if f is f_title else int(size * 0.62)

    for (name, path), (R, C) in zip(items, free):
        page, x, y = slot_xy(R, C)
        im = Image.open(path).convert("RGB").resize((cpx, cpx), Image.LANCZOS)
        page.paste(im, (x, y))
        d = ImageDraw.Draw(page)
        caption(d, name, x, y + cpx + px(0.5), cpx, f_cap)

    if credit:
        R, C = credit[0]
        page, x, y = slot_xy(R, C)
        d = ImageDraw.Draw(page)
        d.text((x, y + px(4)), "Images: Chang et al. 2026, JKAS", fill=(80, 80, 80), font=f_sub)
        d.text((x, y + px(10.5)), f"seowony.github.io/{atlas}.html", fill=(80, 80, 80), font=f_sub)

    # tiny sheet labels in the corner margin so the sheets are easy to arrange
    out = []
    for sy in range(ny):
        for sx in range(nx):
            page = sheets[sy][sx]
            d = ImageDraw.Draw(page)
            lab = f"{sy * nx + sx + 1}/{nx * ny}  row {sy + 1}, column {sx + 1}"
            d.text((W - px(gh / 2) - d.textlength(lab, font=f_tiny), H - px(gv / 2) + px(0.8)), lab,
                   fill=(185, 185, 185), font=f_tiny)
            out.append(page)
    return out, cell_mm, gh, gv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("atlas", choices=list(TITLES))
    ap.add_argument("--out", default="print")
    ap.add_argument("--poster", default="A0", choices=list(PAPER_MM))
    ap.add_argument("--tiles", help="also make a poster over A4 sheets, e.g. 3x2 (columns x rows)")
    ap.add_argument("--poster-jpg", action="store_true", help="also write the poster as a JPEG")
    ap.add_argument("--dpi", type=int, default=300, help="output resolution (pages are drawn at 300 dpi, then resampled)")
    ap.add_argument("--poster-dpi", type=int, help="resolution of the poster (default: --dpi)")
    ap.add_argument("--quality", type=int, default=92, help="JPEG quality inside the PDFs")
    ap.add_argument("--exclude", default="", help="ids to leave out, e.g. ngc2822,ngc5128")
    args = ap.parse_args()
    skip = {x.strip().lower() for x in args.exclude.split(",") if x.strip()}
    items = [it for it in load_items(args.atlas)
             if os.path.basename(it[1])[:-4] not in skip]
    if not items:
        sys.exit(f"No images in assets/{args.atlas}/")
    os.makedirs(args.out, exist_ok=True)

    def out(img, dpi=None):
        dpi = dpi or args.dpi
        if dpi == DPI:
            return img
        f = dpi / DPI
        return img.resize((round(img.width * f), round(img.height * f)), Image.LANCZOS)
    q = args.quality

    if args.tiles:
        nx, ny = (int(v) for v in args.tiles.lower().split("x"))
        sheets, cell, gh, gv = tiles(args.atlas, items, nx, ny)
        tp = os.path.join(args.out, f"{args.atlas}_tiles_{nx}x{ny}.pdf")
        sheets = [out(x) for x in sheets]
        sheets[0].save(tp, save_all=True, append_images=sheets[1:], resolution=args.dpi, quality=q)
        print(f"{tp}: {len(sheets)} A4 sheets ({nx} across, {ny} down), images {cell:.1f} mm, "
              f"gaps {gh:.1f} x {gv:.1f} mm, finished size {210 * nx} x {297 * ny} mm")

    pages = a4_pages(args.atlas, items)
    a4 = os.path.join(args.out, f"{args.atlas}_A4.pdf")
    pages = [out(x) for x in pages]
    pages[0].save(a4, save_all=True, append_images=pages[1:], resolution=args.dpi, quality=q)
    print(f"{a4}: {len(pages)} pages, {len(items)} images")

    pst, cols, rows = poster(args.atlas, items, args.poster)
    pdf = os.path.join(args.out, f"{args.atlas}_poster.pdf")
    pdpi = args.poster_dpi or args.dpi
    pst = out(pst, pdpi)
    pst.save(pdf, resolution=pdpi, quality=q)
    if args.poster_jpg:
        pst.save(os.path.join(args.out, f"{args.atlas}_poster.jpg"), quality=q, dpi=(pdpi, pdpi))
    print(f"{pdf}: {args.poster}, {cols} x {rows} grid, {pst.size[0]} x {pst.size[1]} px")


if __name__ == "__main__":
    main()
