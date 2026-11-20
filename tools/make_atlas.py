#!/usr/bin/env python3
"""
Build image atlases from the KS4 DR1 BVI colour HiPS:

  arp       arp.html       Arp's Atlas of Peculiar Galaxies (Arp 1-338, positions from Sesame)
  ngc       ngc.html       NGC galaxies larger than 2 arcmin (tools/lists/ngc.csv, from OpenNGC)
  phangs    phangs.html    PHANGS-MUSE / PHANGS-JWST Cycle 1 galaxies (tools/lists/phangs.csv)
  clusters  clusters.html  nearby southern galaxy clusters (tools/lists/clusters.csv)

For each object the script
  1. takes the position from the list (or resolves the name through Sesame),
  2. skips it if it is outside the KS4 footprint (HiPS MOC, or a mostly blank cutout),
  3. cuts a TAN image, north up, from the HiPS (CDS hips2fits, or the tiles directly),
  4. writes a large image and a 300 px thumbnail with a scale bar to assets/<collection>/,
  5. fills atlas_template.html and writes <collection>.html.

Run it from the root of the website repository:

    pip install astropy astropy-healpix scipy requests pillow numpy
    python tools/make_atlas.py arp --only 1-20        # quick test
    python tools/make_atlas.py ngc --hips-url <HiPS base URL or local directory>
    python tools/make_atlas.py all --hips-url <HiPS base URL or local directory>

--only takes Arp numbers for arp and row numbers (or ids such as ngc1365) for the lists.
Downloads are cached in assets/<collection>/, so the script can be re-run safely.
Per-object field sizes for Arp can be set in tools/arp_fov.csv ("arp,fov_arcmin").
The NGC and PHANGS lists are rebuilt with tools/make_lists.py.
"""
import argparse, csv, io, json, os, sys, time
from urllib.parse import quote

import numpy as np
import requests
from PIL import Image, ImageDraw

HIPS = "SNU/P/KS4/DR1/colorBVI"
HIPS2FITS = "https://alasky.cds.unistra.fr/hips-image-services/hips2fits"
ALADIN = "https://aladin.cds.unistra.fr/AladinLite/"
THUMB = 300
DEC_LIMIT = -28.0          # only objects south of this declination (deg)
BLANK_LIMIT = 0.35         # skip if more than this fraction of pixels is empty
HOW_TO = ("Click an image to see it larger. The name under each image opens the KS4 <i>BVI</i> colour HiPS in "
          '<a href="https://aladin.cds.unistra.fr/AladinLite/">Aladin Lite</a> so you can pan and zoom around it. ')
KS4_CREDIT = ("Images: KMTNet Synoptic Survey of Southern Sky (KS4) DR1, <i>BVI</i> colour HiPS "
              "(Chang et al. 2026, JKAS).")

COLLECTIONS = {
    "arp": dict(
        title="Arp galaxies in KS4", short="Arp galaxies", big=1024, fov=5.0,
        description="Galaxies from Halton Arp's Atlas of Peculiar Galaxies, seen in the KMTNet Synoptic "
                    "Survey of Southern Sky (KS4) Data Release 1.",
        intro='This gallery shows the galaxies from <a href="https://ned.ipac.caltech.edu/level5/Arp/'
              'Arp_contents.html">Halton Arp\'s Atlas of Peculiar Galaxies</a> that lie south of declination '
              "{dec} and fall inside the footprint of KS4 Data Release 1. " + HOW_TO,
        count="{n} Arp objects south of {dec} are inside the KS4 DR1 footprint.",
        credit=KS4_CREDIT + " Object list: Arp (1966), positions from SIMBAD."),
    "ngc": dict(
        title="NGC galaxies in KS4", short="NGC galaxies", big=800, list="tools/lists/ngc.csv",
        description="Galaxies from the New General Catalogue, seen in the KMTNet Synoptic Survey of "
                    "Southern Sky (KS4) Data Release 1.",
        intro="This gallery shows the galaxies from the New General Catalogue that are larger than two "
              "arcminutes, lie south of declination {dec} and fall inside the footprint of KS4 Data Release 1. "
              + HOW_TO,
        count="{n} of the {total} NGC galaxies south of {dec} are inside the KS4 DR1 footprint.",
        credit=KS4_CREDIT + ' Object list, positions and sizes: <a href="https://github.com/mattiaverga/'
               'OpenNGC">OpenNGC</a> (CC BY-SA 4.0).'),
    "phangs": dict(
        title="PHANGS galaxies in KS4", short="PHANGS galaxies", big=1024, list="tools/lists/phangs.csv",
        description="Nearby star-forming galaxies from the PHANGS surveys, seen in the KMTNet Synoptic "
                    "Survey of Southern Sky (KS4) Data Release 1.",
        intro='This gallery shows the galaxies of the <a href="https://sites.google.com/view/phangs/home">'
              "PHANGS</a>-MUSE and PHANGS-JWST Cycle 1 samples that lie south of declination {dec} and fall "
              "inside the footprint of KS4 Data Release 1. " + HOW_TO,
        count="{n} of the {total} PHANGS galaxies south of {dec} are inside the KS4 DR1 footprint.",
        credit=KS4_CREDIT + ' Positions and sizes: <a href="https://github.com/mattiaverga/OpenNGC">'
               "OpenNGC</a> (CC BY-SA 4.0)."),
    "clusters": dict(
        title="Nearby galaxy clusters in KS4", short="Galaxy clusters", big=1600,
        list="tools/lists/clusters.csv",
        description="Nearby southern galaxy clusters seen in the KMTNet Synoptic Survey of Southern Sky "
                    "(KS4) Data Release 1.",
        intro="This gallery shows nearby galaxy clusters south of declination {dec} that fall inside the "
              "footprint of KS4 Data Release 1. Each image is centred on the brightest cluster galaxy or on the "
              "cluster position in SIMBAD. " + HOW_TO,
        count="{n} of the {total} clusters south of {dec} are inside the KS4 DR1 footprint.",
        credit=KS4_CREDIT + " Positions from SIMBAD and OpenNGC."),
}


def parse_only(text, n_max):
    if not text:
        return list(range(1, n_max + 1))
    nums = []
    for part in text.split(","):
        if part.strip().isdigit():
            nums.append(int(part))
        elif "-" in part and all(x.strip().isdigit() for x in part.split("-")):
            a, b = part.split("-")
            nums += range(int(a), int(b) + 1)
        else:
            nums.append(part.strip().lower())
    return nums


def load_objects(name, cfg, only, dec_limit):
    """List of dicts: id, name, query, ra, dec, fov (ra/dec may be None)."""
    if name == "arp":
        fov = {}
        if os.path.exists("tools/arp_fov.csv"):
            with open("tools/arp_fov.csv") as f:
                for row in csv.DictReader(f):
                    fov[int(row["arp"])] = float(row["fov_arcmin"])
        all_objs = [dict(id=f"arp{n:03d}", name=f"Arp {n}", query=f"Arp {n}", alt=f"ARP {n:03d}",
                         ra=None, dec=None, fov=fov.get(n, cfg["fov"])) for n in range(1, 339)]
    else:
        all_objs = []
        with open(cfg["list"]) as f:
            for row in csv.DictReader(f):
                all_objs.append(dict(id=row["id"], name=row["name"], query=row.get("query") or row["name"],
                                     ra=float(row["ra"]) if row.get("ra") else None,
                                     dec=float(row["dec"]) if row.get("dec") else None,
                                     fov=float(row["fov_arcmin"])))
        n_all = len(all_objs)
        all_objs = [o for o in all_objs if o["dec"] is None or o["dec"] < dec_limit]
        if len(all_objs) < n_all:
            print(f"{name}: {n_all - len(all_objs)} objects north of {dec_limit:g} deg left out")
    sel = parse_only(only, len(all_objs))
    by_id = {o["id"]: o for o in all_objs}
    out = []
    for k in sel:
        if isinstance(k, int) and 1 <= k <= len(all_objs):
            out.append(all_objs[k - 1])
        elif k in by_id:
            out.append(by_id[k])
        else:
            print(f"--only: {k} not in the {name} list")
    return all_objs, out


def resolve(obj, cache):
    key = obj["id"]
    if obj["ra"] is not None:
        return [obj["ra"], obj["dec"]]
    if key in cache:
        return cache[key]
    from astropy.coordinates import SkyCoord
    for name in filter(None, (obj["query"], obj.get("alt"))):
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
TILE_STATS = {"found": 0, "missing": 0, "used": 0}   # used: real tiles placed in cutouts
_HAS = {}


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
            _HAS[key] = False
            arr = np.zeros((tile_size, tile_size, 3), dtype=np.uint8)
        else:
            TILE_STATS["found"] += 1
            _HAS[key] = True
            # JPEG/PNG tiles are stored top row first, the HPX array bottom row first.
            arr = np.asarray(Image.open(io.BytesIO(data)).convert("RGB"))[::-1]
        if len(_TILES) > 400:
            _TILES.clear()
        _TILES[key] = arr
    if _HAS.get(key):
        TILE_STATS["used"] += 1
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


def bar_arcmin(fov_arcmin):
    return 5 if fov_arcmin >= 25 else 1


def scale_bar(img, fov_arcmin):
    w, h = img.size
    px = w / fov_arcmin * bar_arcmin(fov_arcmin)
    m = int(round(w * 0.04))
    t = max(2, int(round(w * 0.006)))
    d = ImageDraw.Draw(img)
    d.rectangle([w - m - px, h - m - t, w - m, h - m], fill=(255, 255, 255))
    return img


def build(name, args):
    cfg = COLLECTIONS[name]
    out_dir = f"assets/{name}"
    page_name = f"{name}.html"
    big_px = cfg["big"]
    os.makedirs(out_dir, exist_ok=True)
    all_objs, todo = load_objects(name, cfg, args.only, args.dec_limit)

    cache_file = os.path.join(out_dir, "positions.json")
    pos = json.load(open(cache_file)) if os.path.exists(cache_file) else {}
    status_file = os.path.join(out_dir, "status.json")
    status = json.load(open(status_file)) if os.path.exists(status_file) else {}
    if name == "arp":   # caches from the first version were keyed "1", "2", ...
        for d in (pos, status):
            for k in [k for k in d if k.isdigit()]:
                d.setdefault(f"arp{int(k):03d}", d.pop(k))

    print(f"== {name}: {len(todo)} objects")
    for obj in todo:
        oid, label = obj["id"], obj["name"]
        big = os.path.join(out_dir, f"{oid}.jpg")
        if os.path.exists(big) and status.get(oid) == "ok":
            continue
        rd = resolve(obj, pos)
        json.dump(pos, open(cache_file, "w"), indent=0)
        if rd is None:
            print(f"{label}: name not resolved")
            status[oid] = "unresolved"
            continue
        obj["ra"], obj["dec"] = rd
        if rd[1] >= args.dec_limit:
            status[oid] = "north"
            continue
        f_arc = obj["fov"]
        if args.hips_url:
            base = args.hips_url.rstrip("/")
            if not args.keep_blank and in_moc(base, rd[0], rd[1]) is False:
                print(f"{label}: outside KS4 (MOC, {rd[0]:.2f} {rd[1]:+.2f})")
                status[oid] = "outside"
                continue
            before = dict(TILE_STATS)
            img = hips_cutout(base, rd[0], rd[1], f_arc, big_px)
            if img is not None and TILE_STATS["used"] == before["used"]:
                print(f"{label}: no tiles found ({rd[0]:.2f} {rd[1]:+.2f}, "
                      f"{TILE_STATS['missing'] - before['missing']} missing)")
                status[oid] = "outside"
                continue
        else:
            img = cutout(rd[0], rd[1], f_arc, big_px)
            time.sleep(args.pause)
        if img is None:
            print(f"{label}: download failed")
            status[oid] = "failed"
            continue
        bf = blank_fraction(img)
        if bf > BLANK_LIMIT and not args.keep_blank:
            print(f"{label}: outside KS4 ({bf:.0%} blank, {rd[0]:.2f} {rd[1]:+.2f})")
            status[oid] = "outside"
            continue
        scale_bar(img, f_arc).save(big, quality=86)
        th = img.resize((THUMB, THUMB), Image.LANCZOS)
        scale_bar(th, f_arc).save(os.path.join(out_dir, f"{oid}_t.jpg"), quality=85)
        status[oid] = "ok"
        print(f"{label}: ok ({rd[0]:.4f}, {rd[1]:+.4f}, {f_arc}')")
        json.dump(status, open(status_file, "w"), indent=0)

    json.dump(status, open(status_file, "w"), indent=0)

    # ---- write the page, in list order ----
    for o in all_objs:
        if o["ra"] is None and o["id"] in pos and pos[o["id"]]:
            o["ra"], o["dec"] = pos[o["id"]]
    ok = [o for o in all_objs if status.get(o["id"]) == "ok" and o["dec"] is not None
          and o["dec"] < args.dec_limit and os.path.exists(os.path.join(out_dir, f"{o['id']}.jpg"))]
    north = sum(1 for o in all_objs if status.get(o["id"]) == "north")
    if north:
        print(f"north of {args.dec_limit:g} deg: {north} (left out)")
    if args.hips_url:
        print(f"Tiles read: {TILE_STATS['found']}, missing: {TILE_STATS['missing']}")
    for k in ("outside", "unresolved", "failed"):
        bad = [o["name"] for o in all_objs if status.get(o["id"]) == k]
        if bad:
            print(f"{k}: {len(bad)} -> {', '.join(bad)}")
    if not ok:
        print(f"No images: {page_name} was not written.")
        return False
    figs = []
    for o in ok:
        ra, dec = o["ra"], o["dec"]
        if ra is None:
            ra, dec = pos[o["id"]]
        f_arc = o["fov"]
        target = quote(f"{ra:.5f} {dec:+.5f}")
        link = f"{ALADIN}?target={target}&fov={max(f_arc * 2, 6) / 60:.3f}&survey={quote(HIPS, safe='')}"
        bar = bar_arcmin(f_arc)
        alt = f"{o['name']} in KS4 BVI colour, {f_arc:g} arcmin field, scale bar {bar} arcmin"
        figs.append(
            f'      <figure>\n'
            f'        <a class="thumb" href="{out_dir}/{o["id"]}.jpg"><img src="{out_dir}/{o["id"]}_t.jpg" '
            f'width="{THUMB}" height="{THUMB}" loading="lazy" alt="{alt}"></a>\n'
            f'        <figcaption><a href="{link}">{o["name"]}</a></figcaption>\n'
            f'      </figure>')
    grid = '    <div class="arp-grid">\n' + "\n".join(figs) + '\n    </div>'
    bars = sorted({bar_arcmin(o["fov"]) for o in ok})
    bar_text = " and ".join({1: "one arcminute", 5: "five arcminutes"}[b] for b in bars)
    dec_txt = f"\u2212{abs(args.dec_limit):g}\u00b0" if args.dec_limit < 0 else f"{args.dec_limit:g}\u00b0"
    intro = cfg["intro"].replace("{dec}", dec_txt) + f"The white bar in the lower-right corner of each image is {bar_text}."
    page = open("atlas_template.html").read()
    for key, val in (("TITLE", cfg["title"]), ("SHORT", cfg["short"]), ("DESCRIPTION", cfg["description"]),
                     ("PAGE", page_name), ("INTRO", intro), ("CREDIT", cfg["credit"]),
                     ("COUNT", cfg["count"].format(n=len(ok), total=len(all_objs), dec=dec_txt)), ("GRID", grid)):
        page = page.replace(f"<!-- {key} -->", val)
    open(page_name, "w").write(page)
    print(f"Wrote {page_name} with {len(ok)} objects.")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("collection", nargs="?", default="arp", choices=list(COLLECTIONS) + ["all"])
    ap.add_argument("--only", help="arp: 1-50 or 81,82,244; lists: row numbers or ids such as ngc1365")
    ap.add_argument("--pause", type=float, default=0.5, help="seconds between hips2fits requests")
    ap.add_argument("--keep-blank", action="store_true", help="keep objects outside the footprint")
    ap.add_argument("--dec-limit", type=float, default=DEC_LIMIT,
                    help=f"keep only objects south of this declination (default {DEC_LIMIT:g})")
    ap.add_argument("--hips-url", help="read tiles from this HiPS (URL or directory) instead of CDS hips2fits")
    args = ap.parse_args()

    if not os.path.exists("atlas_template.html"):
        sys.exit("Run this from the website root (atlas_template.html not found).")
    names = list(COLLECTIONS) if args.collection == "all" else [args.collection]
    if args.only and len(names) > 1:
        sys.exit("--only works with one collection at a time.")
    results = [build(n, args) for n in names]
    if not any(results):
        sys.exit(1)


if __name__ == "__main__":
    main()
