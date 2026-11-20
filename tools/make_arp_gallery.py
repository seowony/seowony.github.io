#!/usr/bin/env python3
"""
Build arp.html: Arp's Atlas of Peculiar Galaxies as seen in KS4 DR1.

For each of the 338 Arp objects this script
  1. resolves the position by name through Sesame (SIMBAD/NED),
  2. asks CDS hips2fits for a BVI colour cutout from the KS4 DR1 HiPS,
  3. skips objects that fall outside the KS4 footprint (mostly blank cutouts),
  4. writes a 1024 px image and a 300 px thumbnail, each with a 1-arcmin bar,
  5. fills arp_template.html and writes arp.html.

Run it from the root of the website repository:

    pip install astropy astropy-healpix scipy requests pillow numpy
    python tools/make_arp_gallery.py --only 1-20      # quick test
    python tools/make_arp_gallery.py                  # all 338 objects

By default the cutouts come from the public KS4 HiPS at CDS (hips2fits).
To read the tiles straight from another copy of the HiPS, for example the
internal SNU server, give its base URL (or a local directory):

    python tools/make_arp_gallery.py \
        --hips-url <HiPS base URL or local directory>

Downloads are cached in assets/arp/, so the script can be re-run safely.
Per-object field sizes can be set in tools/arp_fov.csv ("arp,fov_arcmin").
"""
import argparse, csv, io, json, os, sys, time
from urllib.parse import quote

import numpy as np
import requests
from PIL import Image, ImageDraw

HIPS = "SNU/P/KS4/DR1/colorBVI"
HIPS2FITS = "https://alasky.cds.unistra.fr/hips-image-services/hips2fits"
ALADIN = "https://aladin.cds.unistra.fr/AladinLite/"
OUT = "assets/arp"
BIG, THUMB = 1024, 300
DEFAULT_FOV = 5.0          # arcmin
BLANK_LIMIT = 0.35         # skip if more than this fraction of pixels is empty


def parse_only(text):
    if not text:
        return list(range(1, 339))
    nums = []
    for part in text.split(","):
        if "-" in part:
            a, b = part.split("-")
            nums += range(int(a), int(b) + 1)
        else:
            nums.append(int(part))
    return nums


def resolve(n, cache):
    key = str(n)
    if key in cache:
        return cache[key]
    from astropy.coordinates import SkyCoord
    for name in (f"Arp {n}", f"ARP {n:03d}"):
        try:
            c = SkyCoord.from_name(name)
            cache[key] = [c.ra.deg, c.dec.deg]
            return cache[key]
        except Exception:
            continue
    cache[key] = None
    return None


def cutout(ra, dec, fov_arcmin, size):
    params = dict(hips=HIPS, width=size, height=size, fov=fov_arcmin / 60.0,
                  projection="TAN", coordsys="icrs", ra=ra, dec=dec, format="jpg")
    for attempt in range(4):
        try:
            r = requests.get(HIPS2FITS, params=params, timeout=90)
            if r.status_code == 200 and r.headers.get("content-type", "").startswith("image"):
                return Image.open(io.BytesIO(r.content)).convert("RGB")
        except requests.RequestException:
            pass
        time.sleep(3 * (attempt + 1))
    return None


# ---------------------------------------------------------------------------
# Direct HiPS reader (used with --hips-url)
# ---------------------------------------------------------------------------
_PROPS = {}
_TILES = {}
_MOC = {}
TILE_STATS = {"found": 0, "missing": 0}


def _fetch(url_or_path, binary=True):
    if url_or_path.startswith(("http://", "https://")):
        r = requests.get(url_or_path, timeout=60)
        if r.status_code != 200:
            return None
        return r.content if binary else r.text
    if not os.path.exists(url_or_path):
        return None
    with open(url_or_path, "rb" if binary else "r") as f:
        return f.read()


def hips_properties(base):
    if base not in _PROPS:
        text = _fetch(base + "/properties", binary=False)
        if text is None:
            sys.exit(f"Cannot read {base}/properties")
        props = {}
        for line in text.splitlines():
            if "=" in line and not line.lstrip().startswith("#"):
                k, v = line.split("=", 1)
                props[k.strip()] = v.strip()
        _PROPS[base] = props
    return _PROPS[base]


def in_moc(base, ra, dec):
    """True/False from the HiPS Moc.fits, or None if there is no MOC to check against."""
    if base not in _MOC:
        _MOC[base] = None
        data = _fetch(base + "/Moc.fits")
        if data is not None:
            from astropy.io import fits
            with fits.open(io.BytesIO(data)) as h:
                uniq = np.asarray(h[1].data.field(0), dtype=np.int64)
            order = (np.floor(np.log2(uniq) / 2) - 1).astype(int)
            _MOC[base] = (set(uniq.tolist()), sorted(set(order.tolist())))
    if _MOC[base] is None:
        return None
    from astropy_healpix import HEALPix
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    uniq, orders = _MOC[base]
    frame = "galactic" if hips_properties(base).get("hips_frame", "equatorial").startswith("gal") else "icrs"
    c = SkyCoord(ra * u.deg, dec * u.deg)
    for o in orders:
        ipix = int(HEALPix(nside=2 ** o, order="nested", frame=frame).skycoord_to_healpix(c))
        if 4 * 4 ** o + ipix in uniq:
            return True
    return False


def _hpx_wcs(level, frame, tile_size):
    """WCS of the whole HiPS level as one HPX image (same layout as the reproject package)."""
    from astropy.wcs import WCS
    nside = 2 ** level
    size = 5 * nside * tile_size
    res = 45 / tile_size / 2 ** level
    w = WCS(naxis=2)
    w.wcs.ctype = ["GLON-HPX", "GLAT-HPX"] if frame == "galactic" else ["RA---HPX", "DEC--HPX"]
    w.wcs.crval = [0.0, 0.0]
    w.wcs.cd = [[-res, -res], [res, -res]]
    w.wcs.set_pv([(2, 1, 4), (2, 2, 3)])
    w.wcs.crpix = [size / 2 + 0.5, size / 2 + 0.5]
    w.pixel_shape = (size, size)
    return w


def _tile(base, level, index, ext, tile_size):
    key = (base, level, index)
    if key not in _TILES:
        path = f"{base}/Norder{level}/Dir{(index // 10000) * 10000}/Npix{index}.{ext}"
        data = _fetch(path)
        if data is None:
            TILE_STATS["missing"] += 1
            arr = np.zeros((tile_size, tile_size, 3), dtype=np.uint8)
        else:
            TILE_STATS["found"] += 1
            # JPEG/PNG tiles are stored top row first, the HPX array bottom row first.
            arr = np.asarray(Image.open(io.BytesIO(data)).convert("RGB"))[::-1]
        if len(_TILES) > 400:
            _TILES.clear()
        _TILES[key] = arr
    return _TILES[key]


def hips_cutout(base, ra, dec, fov_arcmin, size):
    """TAN cutout centred on (ra, dec) built from the HiPS tiles."""
    import math
    from astropy.coordinates import SkyCoord
    from astropy.wcs import WCS
    from astropy_healpix import HEALPix
    from scipy.ndimage import map_coordinates
    import astropy.units as u

    p = hips_properties(base)
    tw = int(p.get("hips_tile_width", 512))
    max_order = int(p["hips_order"])
    frame = "galactic" if p.get("hips_frame", "equatorial").startswith("gal") else "icrs"
    fmts = p.get("hips_tile_format", "jpeg").split()
    ext = "jpg" if any(f in ("jpeg", "jpg") for f in fmts) else "png"

    out_res = fov_arcmin / 60.0 / size
    level = min(max_order, max(0, math.ceil(math.log2(45 / (tw * out_res)))))
    hpx = _hpx_wcs(level, frame, tw)
    hp = HEALPix(nside=2 ** level, order="nested", frame=frame)

    out = WCS(naxis=2)
    out.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    out.wcs.crval = [ra, dec]
    out.wcs.crpix = [size / 2 + 0.5, size / 2 + 0.5]
    out.wcs.cdelt = [-out_res, out_res]
    yy, xx = np.mgrid[0:size, 0:size]
    world = out.pixel_to_world(xx, yy)
    if frame == "galactic":
        world = world.galactic
    x, y = hpx.world_to_pixel(world)
    good = np.isfinite(x) & np.isfinite(y)
    if not good.any():
        return None
    tx = np.floor((x[good] + 0.5) / tw).astype(int)
    ty = np.floor((y[good] + 0.5) / tw).astype(int)
    tx0, tx1, ty0, ty1 = tx.min(), tx.max(), ty.min(), ty.max()
    mosaic = np.zeros(((ty1 - ty0 + 1) * tw, (tx1 - tx0 + 1) * tw, 3), dtype=np.uint8)
    for jx in range(tx0, tx1 + 1):
        for jy in range(ty0, ty1 + 1):
            px = np.array([jx * tw + 0.25 * tw, jx * tw + 0.75 * tw])
            py = np.array([jy * tw + 0.25 * tw, jy * tw + 0.75 * tw])
            c = hpx.pixel_to_world(px, py)
            try:
                idx = hp.skycoord_to_healpix(c)
            except Exception:
                continue
            idx = idx[idx >= 0]
            if idx.size == 0:
                continue
            arr = _tile(base, level, int(idx.max()), ext, tw)
            mosaic[(jy - ty0) * tw:(jy - ty0 + 1) * tw, (jx - tx0) * tw:(jx - tx0 + 1) * tw] = arr[:tw, :tw]
    my = np.where(good, y - ty0 * tw, -1)
    mx = np.where(good, x - tx0 * tw, -1)
    img = np.zeros((size, size, 3), dtype=np.uint8)
    for k in range(3):
        img[..., k] = np.clip(map_coordinates(mosaic[..., k].astype(float), [my, mx],
                                              order=1, mode="constant", cval=0), 0, 255)
    # output row 0 is the bottom of the sky image; flip so north is up
    return Image.fromarray(img[::-1])


def blank_fraction(img):
    a = np.asarray(img).astype(int).sum(axis=2)
    return float((a <= 6).mean())


def scale_bar(img, fov_arcmin):
    w, h = img.size
    px = w / fov_arcmin            # pixels per arcmin
    m = int(round(w * 0.04))
    t = max(2, int(round(w * 0.006)))
    d = ImageDraw.Draw(img)
    d.rectangle([w - m - px, h - m - t, w - m, h - m], fill=(255, 255, 255))
    return img


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="e.g. 1-50 or 81,82,244")
    ap.add_argument("--pause", type=float, default=0.5, help="seconds between requests")
    ap.add_argument("--keep-blank", action="store_true", help="keep objects outside the footprint")
    ap.add_argument("--hips-url", help="read tiles from this HiPS (URL or directory) instead of CDS hips2fits")
    args = ap.parse_args()

    if not os.path.exists("arp_template.html"):
        sys.exit("Run this from the website root (arp_template.html not found).")
    os.makedirs(OUT, exist_ok=True)

    fov = {}
    if os.path.exists("tools/arp_fov.csv"):
        with open("tools/arp_fov.csv") as f:
            for row in csv.DictReader(f):
                fov[int(row["arp"])] = float(row["fov_arcmin"])

    cache_file = os.path.join(OUT, "positions.json")
    pos = json.load(open(cache_file)) if os.path.exists(cache_file) else {}
    status_file = os.path.join(OUT, "status.json")
    status = json.load(open(status_file)) if os.path.exists(status_file) else {}

    for n in parse_only(args.only):
        big = os.path.join(OUT, f"arp{n:03d}.jpg")
        if os.path.exists(big) and status.get(str(n)) == "ok":
            continue
        rd = resolve(n, pos)
        json.dump(pos, open(cache_file, "w"), indent=0)
        if rd is None:
            print(f"Arp {n}: name not resolved")
            status[str(n)] = "unresolved"
            continue
        f_arc = fov.get(n, DEFAULT_FOV)
        if args.hips_url:
            base = args.hips_url.rstrip("/")
            if not args.keep_blank and in_moc(base, rd[0], rd[1]) is False:
                print(f"Arp {n}: outside KS4 (MOC, {rd[0]:.2f} {rd[1]:+.2f})")
                status[str(n)] = "outside"
                continue
            before = dict(TILE_STATS)
            img = hips_cutout(base, rd[0], rd[1], f_arc, BIG)
            if img is not None and TILE_STATS["found"] == before["found"]:
                print(f"Arp {n}: no tiles found ({rd[0]:.2f} {rd[1]:+.2f}, "
                      f"{TILE_STATS['missing'] - before['missing']} missing)")
                status[str(n)] = "outside"
                continue
        else:
            img = cutout(rd[0], rd[1], f_arc, BIG)
            time.sleep(args.pause)
        if img is None:
            print(f"Arp {n}: download failed")
            status[str(n)] = "failed"
            continue
        bf = blank_fraction(img)
        if bf > BLANK_LIMIT and not args.keep_blank:
            print(f"Arp {n}: outside KS4 ({bf:.0%} blank, {rd[0]:.2f} {rd[1]:+.2f})")
            status[str(n)] = "outside"
            continue
        scale_bar(img, f_arc).save(big, quality=88)
        th = img.resize((THUMB, THUMB), Image.LANCZOS)
        scale_bar(th, f_arc).save(os.path.join(OUT, f"arp{n:03d}_t.jpg"), quality=85)
        status[str(n)] = "ok"
        print(f"Arp {n}: ok ({rd[0]:.4f}, {rd[1]:+.4f}, {f_arc}')")
        json.dump(status, open(status_file, "w"), indent=0)

    json.dump(status, open(status_file, "w"), indent=0)

    # ---- write arp.html ----
    ok = sorted(int(k) for k, v in status.items() if v == "ok"
                and os.path.exists(os.path.join(OUT, f"arp{int(k):03d}.jpg")))
    if args.hips_url:
        print(f"\nTiles read: {TILE_STATS['found']}, missing: {TILE_STATS['missing']}")
    if not ok:
        for k in ("outside", "unresolved", "failed"):
            bad = sorted(int(x) for x, v in status.items() if v == k)
            if bad:
                print(f"{k}: {len(bad)} -> {bad}")
        sys.exit("\nNo images: arp.html was not written.")
    figs = []
    for n in ok:
        ra, dec = pos[str(n)]
        f_arc = fov.get(n, DEFAULT_FOV)
        target = quote(f"{ra:.5f} {dec:+.5f}")
        link = f"{ALADIN}?target={target}&fov={max(f_arc * 2, 6) / 60:.3f}&survey={quote(HIPS, safe='')}"
        figs.append(
            f'      <figure>\n'
            f'        <a class="thumb" href="{OUT}/arp{n:03d}.jpg"><img src="{OUT}/arp{n:03d}_t.jpg" '
            f'width="{THUMB}" height="{THUMB}" loading="lazy" alt="Arp {n} in KS4 BVI colour"></a>\n'
            f'        <figcaption><a href="{link}">Arp {n}</a></figcaption>\n'
            f'      </figure>')
    grid = ('    <div class="arp-grid">\n' + "\n".join(figs) + '\n    </div>') if figs else \
           '    <div class="arp-empty">No images yet.</div>'
    count = f"{len(ok)} of the 338 Arp objects are inside the KS4 DR1 footprint."
    page = open("arp_template.html").read().replace("<!-- ARP_GRID -->", grid) \
                                         .replace("<!-- ARP_COUNT -->", count)
    open("arp.html", "w").write(page)
    print(f"\nWrote arp.html with {len(ok)} objects.")
    for k in ("outside", "unresolved", "failed"):
        bad = sorted(int(x) for x, v in status.items() if v == k)
        if bad:
            print(f"{k}: {len(bad)} -> {bad}")


if __name__ == "__main__":
    main()
