# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "pymol-open-source-whl",
#     "pillow>=10.0",
#     "numpy>=1.24",
# ]
# ///
"""
03_render_pymol_3d_poses.py
Renders publication-grade, 3D ray-traced PyMOL structural binding-mode panels (300 DPI)
and saves reusable .pse sessions for Protocol B:
  - Panel A: DPP-4 (4A5S) crystal N7F-901 vs. redocked N7F-901 (RMSD = 0.780 Å)
  - Panel B: DPP-4 (4A5S) S1/S2/S2' pocket overlay (Teneligliptin, Naringenin, Berberine)
  - Panel C: CD36 (5LGD) ectodomain hydrophobic tunnel cutaway (LYS164 mouth vs. PLM-512,
             Naringenin upper chamber, Linagliptin & Teneligliptin mid-tunnel)
  - Panel D: CD36 (5LGD) zoomed pocket contacts flanking shared pivot GLU335
"""
import os
import sys
import glob
from PIL import Image, ImageDraw, ImageFont

# Headless setup: use OSMesa on Linux; Windows uses native offscreen PyMOL ray tracer
if os.name != "nt":
    os.environ["PYOPENGL_PLATFORM"] = "osmesa"

import pymol  # pytype: disable=import-error
pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()

from pymol import cmd, util  # pytype: disable=import-error

PROT_B_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(PROT_B_DIR)
STRUCT_DIR = os.path.join(ROOT_DIR, "structures")
POSES_DPP4 = os.path.join(PROT_B_DIR, "results", "poses", "DPP4_4A5S")
POSES_CD36 = os.path.join(PROT_B_DIR, "results", "poses", "CD36_5LGD")
OUT_DIR = os.path.join(PROT_B_DIR, "figures", "pymol_3d")
os.makedirs(OUT_DIR, exist_ok=True)


def extract_model1_pdb(pdbqt_path, out_pdb, resname="LIG"):
    """Converts Model 1 of a Vina PDBQT pose to clean standard PDB so PyMOL infers covalent bonds cleanly."""
    lines = []
    in_m1 = False
    has_model = False
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
        raw = f.readlines()
    for l in raw:
        if l.startswith("MODEL"):
            has_model = True
            if l.strip().split()[-1] == "1":
                in_m1 = True
            else:
                break
        elif l.startswith("ENDMDL") and in_m1:
            break
        elif (in_m1 or not has_model) and l.startswith(("ATOM  ", "HETATM")):
            aname = l[12:16].strip()
            adtype = l[77:79].strip() if len(l) >= 78 else ""
            if adtype in ("H", "HD") or aname.startswith("H"):
                continue
            el_map = {"A": "C", "C": "C", "N": "N", "NA": "N", "O": "O", "OA": "O", "S": "S", "SA": "S", "F": "F", "CL": "CL", "BR": "BR", "P": "P"}
            el = el_map.get(adtype.upper(), "".join(c for c in aname if c.isalpha())[:2].upper() or "C")
            # Build clean HETATM line with standard right-justified element column (77-78)
            pdb_line = f"HETATM{l[6:12]}{l[12:16]} {resname:>3s} L   1    {l[30:54]}  1.00 20.00          {el:>2s}\n"
            lines.append(pdb_line)
    with open(out_pdb, "w", encoding="utf-8") as fo:
        fo.writelines(lines)
        fo.write("END\n")
    return out_pdb


def find_pose(pdir, pat):
    m = glob.glob(os.path.join(pdir, pat))
    if not m:
        raise FileNotFoundError(f"Pose not found: {os.path.join(pdir, pat)}")
    return m[0]


def setup_pub_style():
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("orthoscopic", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("ambient", 0.35)
    cmd.set("direct", 0.65)
    cmd.set("specular", 0.25)
    cmd.set("depth_cue", 0)
    cmd.set("fog", 0)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("stick_radius", 0.18)
    cmd.set("dash_radius", 0.07)
    cmd.set("dash_gap", 0.30)
    cmd.set("dash_length", 0.22)
    cmd.set("two_sided_lighting", 1)
    cmd.set("label_font_id", 7)
    cmd.set("label_size", 15)
    cmd.set("label_color", "black")
    cmd.set("label_outline_color", "white")
    cmd.set("label_bg_color", "white")
    cmd.set("label_bg_transparency", 0.15)
    cmd.set("label_shadow_mode", 0)
    cmd.set("label_position", [0.0, 0.8, 4.0])


def render_panel_a_dpp4_redock(tmp_dir):
    cmd.reinitialize()
    setup_pub_style()
    rec_p = os.path.join(STRUCT_DIR, "4A5S_DPP4_clean.pdb")
    crys_p = os.path.join(STRUCT_DIR, "N7F901_crystal.pdb")
    red_p = extract_model1_pdb(os.path.join(POSES_DPP4, "N7F901_redocked_4A5S.pdbqt"), os.path.join(tmp_dir, "n7f_redock.pdb"), "RED")

    cmd.load(rec_p, "dpp4")
    cmd.load(crys_p, "n7f_crys")
    cmd.load(red_p, "n7f_red")
    assert cmd.count_atoms("all") > 0, "Failed to load Panel A structures!"

    cmd.hide("everything", "all")
    cmd.show("cartoon", "dpp4")
    cmd.color("gray85", "dpp4")
    cmd.set("cartoon_transparency", 0.35, "dpp4")

    cmd.show("surface", "dpp4")
    cmd.set("surface_color", "white", "dpp4")
    cmd.set("transparency", 0.78, "dpp4")

    cmd.select("pocket_res", "dpp4 and resi 125+205+206+357+547+630+631+662")
    cmd.show("sticks", "pocket_res")
    cmd.set("stick_radius", 0.13, "pocket_res")
    cmd.color("gray65", "pocket_res and elem C")
    util.cnc("pocket_res")
    cmd.label("dpp4 and resi 205+206+630+631+662+357+547 and name CA", "'%s%s' % (resn.capitalize(), resi)")

    cmd.show("sticks", "n7f_crys")
    cmd.color("cyan", "n7f_crys and elem C")
    util.cnc("n7f_crys")
    cmd.set("stick_radius", 0.20, "n7f_crys")

    cmd.show("sticks", "n7f_red")
    cmd.color("marine", "n7f_red and elem C")
    util.cnc("n7f_red")
    cmd.set("stick_radius", 0.16, "n7f_red")

    cmd.distance("hb_n7f", "n7f_red", "dpp4 and resi 205+206+631+662", cutoff=3.5, mode=2)
    cmd.color("firebrick", "hb_n7f")
    cmd.hide("labels", "hb_n7f")

    cmd.orient("n7f_crys | pocket_res")
    cmd.zoom("n7f_crys | pocket_res", buffer=3.2)
    out_png = os.path.join(OUT_DIR, "panelA_dpp4_n7f_redock.png")
    cmd.png(out_png, width=1400, height=1050, dpi=300, ray=1)
    cmd.save(os.path.join(OUT_DIR, "panelA_dpp4_n7f_redock.pse"))
    print(f"[PyMOL Panel A] Saved {out_png} ({cmd.count_atoms('all')} atoms)")
    return out_png


def render_panel_b_dpp4_leads(tmp_dir):
    cmd.reinitialize()
    setup_pub_style()
    rec_p = os.path.join(STRUCT_DIR, "4A5S_DPP4_clean.pdb")
    ten_p = extract_model1_pdb(find_pose(POSES_DPP4, "*_Teneligliptin_DPP4.pdbqt"), os.path.join(tmp_dir, "ten_dpp4.pdb"), "TEN")
    nar_p = extract_model1_pdb(find_pose(POSES_DPP4, "*_Naringenin_DPP4.pdbqt"), os.path.join(tmp_dir, "nar_dpp4.pdb"), "NAR")
    ber_p = extract_model1_pdb(find_pose(POSES_DPP4, "*_Berberine_DPP4.pdbqt"), os.path.join(tmp_dir, "ber_dpp4.pdb"), "BER")

    cmd.load(rec_p, "dpp4")
    cmd.load(ten_p, "ten")
    cmd.load(nar_p, "nar")
    cmd.load(ber_p, "ber")
    assert cmd.count_atoms("all") > 0, "Failed to load Panel B structures!"

    cmd.hide("everything", "all")
    cmd.show("cartoon", "dpp4")
    cmd.color("gray85", "dpp4")
    cmd.set("cartoon_transparency", 0.35, "dpp4")
    cmd.show("surface", "dpp4")
    cmd.set("surface_color", "white", "dpp4")
    cmd.set("transparency", 0.78, "dpp4")

    cmd.select("pocket_res", "dpp4 and resi 125+205+206+357+547+630+631+662+752")
    cmd.show("sticks", "pocket_res")
    cmd.set("stick_radius", 0.12, "pocket_res")
    cmd.color("gray65", "pocket_res and elem C")
    util.cnc("pocket_res")
    cmd.label("dpp4 and resi 125+205+206+357+631+662 and name CA", "'%s%s' % (resn.capitalize(), resi)")
    cmd.label("dpp4 and resi 547 and name OH", "'Tyr547'")
    cmd.label("dpp4 and resi 752 and name OH", "'Tyr752'")

    for obj, col in [("ten", "teal"), ("nar", "forest"), ("ber", "brightorange")]:
        cmd.show("sticks", obj)
        cmd.color(col, f"{obj} and elem C")
        util.cnc(obj)
        cmd.set("stick_radius", 0.19, obj)

    cmd.distance("hb_nar", "nar", "dpp4 and resi 206+357+631", cutoff=3.5, mode=2)
    cmd.color("firebrick", "hb_nar")
    cmd.hide("labels", "hb_nar")
    cmd.distance("hb_ten", "ten", "dpp4 and resi 125+547+631+662", cutoff=3.5, mode=2)
    cmd.color("deeppurple", "hb_ten")
    cmd.hide("labels", "hb_ten")

    cmd.orient("ten | nar | ber | pocket_res")
    cmd.zoom("ten | nar | ber | (dpp4 and resi 752 and name OH)", buffer=4.8)
    out_png = os.path.join(OUT_DIR, "panelB_dpp4_tier1_leads.png")
    cmd.png(out_png, width=1400, height=1050, dpi=300, ray=1)
    cmd.save(os.path.join(OUT_DIR, "panelB_dpp4_tier1_leads.pse"))
    print(f"[PyMOL Panel B] Saved {out_png} ({cmd.count_atoms('all')} atoms)")
    return out_png


def render_panel_c_cd36_tunnel(tmp_dir):
    cmd.reinitialize()
    setup_pub_style()
    rec_p = os.path.join(STRUCT_DIR, "5LGD_CD36_clean.pdb")
    plm_p = os.path.join(STRUCT_DIR, "PLM512_crystal.pdb")
    nar_p = extract_model1_pdb(find_pose(POSES_CD36, "*_Naringenin_CD36.pdbqt"), os.path.join(tmp_dir, "nar_cd36.pdb"), "NAR")
    lin_p = extract_model1_pdb(find_pose(POSES_CD36, "*_Linagliptin_CD36.pdbqt"), os.path.join(tmp_dir, "lin_cd36.pdb"), "LIN")
    ten_p = extract_model1_pdb(find_pose(POSES_CD36, "*_Teneligliptin_CD36.pdbqt"), os.path.join(tmp_dir, "ten_cd36.pdb"), "TEN")

    cmd.load(rec_p, "cd36")
    cmd.load(plm_p, "plm512")
    cmd.load(nar_p, "nar")
    cmd.load(lin_p, "lin")
    cmd.load(ten_p, "ten")
    assert cmd.count_atoms("all") > 0, "Failed to load Panel C structures!"

    cmd.hide("everything", "all")
    cmd.show("cartoon", "cd36")
    cmd.color("gray85", "cd36")
    cmd.set("cartoon_transparency", 0.35, "cd36")

    cmd.show("surface", "cd36")
    cmd.set("surface_color", "palecyan", "cd36")
    cmd.set("transparency", 0.72, "cd36")
    cmd.set("surface_cavity_mode", 1)

    # Highlight distal entrance LYS164 (red) and central pivot GLU335 (navy)
    cmd.select("lys164", "cd36 and resi 164")
    cmd.show("sticks", "lys164")
    cmd.color("firebrick", "lys164 and elem C")
    util.cnc("lys164")
    cmd.label("cd36 and resi 164 and name NZ", "'Lys164 (Gate)'")
    cmd.select("glu335", "cd36 and resi 335")
    cmd.show("sticks", "glu335")
    cmd.color("marine", "glu335 and elem C")
    util.cnc("glu335")
    cmd.label("cd36 and resi 335 and name CA", "'Glu335'")

    for obj, col, rad in [("plm512", "skyblue", 0.22), ("nar", "forest", 0.21), ("lin", "teal", 0.19), ("ten", "violetpurple", 0.19)]:
        cmd.show("sticks", obj)
        cmd.color(col, f"{obj} and elem C")
        util.cnc(obj)
        cmd.set("stick_radius", rad, obj)

    cmd.orient("plm512 | nar | lin | ten | lys164")
    cmd.zoom("plm512 | nar | lin | ten | lys164", buffer=4.0)
    out_png = os.path.join(OUT_DIR, "panelC_cd36_tunnel_overview.png")
    cmd.png(out_png, width=1400, height=1050, dpi=300, ray=1)
    cmd.save(os.path.join(OUT_DIR, "panelC_cd36_tunnel_overview.pse"))
    print(f"[PyMOL Panel C] Saved {out_png} ({cmd.count_atoms('all')} atoms)")
    return out_png


def render_panel_d_cd36_zoom(tmp_dir):
    cmd.reinitialize()
    setup_pub_style()
    rec_p = os.path.join(STRUCT_DIR, "5LGD_CD36_clean.pdb")
    plm_p = os.path.join(STRUCT_DIR, "PLM512_crystal.pdb")
    nar_p = os.path.join(tmp_dir, "nar_cd36.pdb")
    lin_p = os.path.join(tmp_dir, "lin_cd36.pdb")
    ten_p = os.path.join(tmp_dir, "ten_cd36.pdb")

    cmd.load(rec_p, "cd36")
    cmd.load(plm_p, "plm512")
    cmd.load(nar_p, "nar")
    cmd.load(lin_p, "lin")
    cmd.load(ten_p, "ten")
    assert cmd.count_atoms("all") > 0, "Failed to load Panel D structures!"

    cmd.hide("everything", "all")
    cmd.show("cartoon", "cd36")
    cmd.color("gray85", "cd36")
    cmd.set("cartoon_transparency", 0.45, "cd36")

    cmd.select("tunnel_res", "cd36 and resi 141+145+195+203+238+269+334+335+337")
    cmd.show("sticks", "tunnel_res")
    cmd.set("stick_radius", 0.13, "tunnel_res")
    cmd.color("gray65", "tunnel_res and elem C")
    util.cnc("tunnel_res")
    cmd.label("cd36 and resi 141+195+238+269+334+335+337 and name CA", "'%s%s' % (resn.capitalize(), resi)")

    for obj, col, rad in [("plm512", "skyblue", 0.18), ("nar", "forest", 0.21), ("lin", "teal", 0.18), ("ten", "violetpurple", 0.18)]:
        cmd.show("sticks", obj)
        cmd.color(col, f"{obj} and elem C")
        util.cnc(obj)
        cmd.set("stick_radius", rad, obj)

    cmd.distance("hb_nar_c", "nar", "cd36 and resi 141+145+195+337", cutoff=3.5, mode=2)
    cmd.distance("hb_ten_c", "ten", "cd36 and resi 238+269", cutoff=3.5, mode=2)
    cmd.distance("hb_lin_c", "lin", "cd36 and resi 269", cutoff=3.5, mode=2)
    for dname in ("hb_nar_c", "hb_ten_c", "hb_lin_c"):
        cmd.color("firebrick", dname)
        cmd.hide("labels", dname)

    cmd.orient("plm512 | nar | lin | ten")
    cmd.zoom("plm512 | nar | lin | ten", buffer=2.6)
    out_png = os.path.join(OUT_DIR, "panelD_cd36_tunnel_contacts.png")
    cmd.png(out_png, width=1400, height=1050, dpi=300, ray=1)
    cmd.save(os.path.join(OUT_DIR, "panelD_cd36_tunnel_contacts.pse"))
    print(f"[PyMOL Panel D] Saved {out_png} ({cmd.count_atoms('all')} atoms)")
    return out_png


def compose_2x2(pA, pB, pC, pD, out_path):
    imgs = [Image.open(p).convert("RGB") for p in (pA, pB, pC, pD)]
    w, h = imgs[0].size
    hdr_h = 90
    leg_h = 112
    pad = 24
    W = w * 2 + pad * 3
    H = (h + hdr_h) * 2 + leg_h + pad * 4
    canvas = Image.new("RGB", (W, H), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)

    try:
        font_title = ImageFont.truetype("arialbd.ttf", 30)
        font_sub = ImageFont.truetype("arial.ttf", 22)
        font_leg = ImageFont.truetype("arialbd.ttf", 24)
    except Exception:
        font_title = ImageFont.load_default()
        font_sub = font_title
        font_leg = font_title

    titles = [
        ("A. DPP-4 (4A5S) Co-Crystal Superposition", "Crystal N7F-901 (cyan) vs. Redocked (marine; Hungarian RMSD = 0.780 Å)"),
        ("B. DPP-4 (4A5S) Tier 1 S1/S2/S2' Pocket Overlay", "Teneligliptin (teal), Naringenin (green), Berberine (orange); H-bonds dashed"),
        ("C. Hepatic CD36 (5LGD) Hydrophobic Tunnel Cutaway", "Lys164 gate (red) vs. PLM-512 (sky-blue), Naringenin (green), Gliptins flanking Glu335 (navy)"),
        ("D. CD36 (5LGD) Adjacent Segment Pocket Contacts", "Upper chamber (Ala141/145, Thr195, Arg337) & mid-tunnel (Ser269, Tyr238, Lys334, Pro203)"),
    ]

    positions = [
        (pad, pad),
        (pad * 2 + w, pad),
        (pad, pad * 2 + h + hdr_h),
        (pad * 2 + w, pad * 2 + h + hdr_h),
    ]
    for (x, y), im, (t1, t2) in zip(positions, imgs, titles):
        draw.rectangle([x, y, x + w, y + hdr_h - 6], fill=(241, 245, 249), outline=(203, 213, 225), width=2)
        draw.text((x + 16, y + 10), t1, fill=(11, 60, 93), font=font_title)
        draw.text((x + 16, y + 48), t2, fill=(71, 85, 105), font=font_sub)
        canvas.paste(im, (x, y + hdr_h))
        draw.rectangle([x, y + hdr_h, x + w, y + hdr_h + h], outline=(203, 213, 225), width=2)

    # Bottom crisp 2-row centered color key banner with measured text widths (0 overlap)
    ly = H - leg_h - pad
    draw.rectangle([pad, ly, W - pad, ly + leg_h], fill=(248, 250, 252), outline=(148, 163, 184), width=2)
    row1 = [
        ((14, 165, 233), "Crystal Controls (N7F-901 / PLM-512)"),
        ((30, 58, 138), "Redocked N7F-901 / Glu335 Pivot"),
        ((15, 118, 110), "Teneligliptin (4A5S) / Linagliptin (5LGD)"),
        ((124, 58, 237), "Teneligliptin (5LGD)"),
    ]
    row2 = [
        ((21, 128, 61), "Naringenin (Tier 1A-NP)"),
        ((217, 119, 6), "Berberine (Tier 1B-Cation)"),
        ((185, 28, 28), "HD-Verified H-Bonds (≤ 3.5 Å) / Lys164 Gate"),
        ((120, 120, 120), "Pocket Side Chains (Element-Colored)"),
    ]
    for r_idx, row_items in enumerate((row1, row2)):
        item_widths = []
        for _, lbl in row_items:
            bb = draw.textbbox((0, 0), lbl, font=font_leg)
            item_widths.append(36 + (bb[2] - bb[0]))
        gap = 54
        total_w = sum(item_widths) + gap * (len(row_items) - 1)
        cx = (W - total_w) // 2
        cy = ly + (16 if r_idx == 0 else 64)
        for (col, lbl), iw in zip(row_items, item_widths):
            draw.rectangle([cx, cy + 3, cx + 24, cy + 27], fill=col, outline=(15, 23, 42), width=1)
            draw.text((cx + 34, cy + 2), lbl, fill=(15, 23, 42), font=font_leg)
            cx += iw + gap

    canvas.save(out_path, dpi=(300, 300))
    print(f"[SUCCESS] Composed 300 DPI PyMOL 3D figure: {out_path}")


def main():
    tmp_dir = os.path.join(OUT_DIR, "_tmp_pdb")
    os.makedirs(tmp_dir, exist_ok=True)
    pA = render_panel_a_dpp4_redock(tmp_dir)
    pB = render_panel_b_dpp4_leads(tmp_dir)
    pC = render_panel_c_cd36_tunnel(tmp_dir)
    pD = render_panel_d_cd36_zoom(tmp_dir)
    out_fig = os.path.join(PROT_B_DIR, "figures", "fig6_pymol_3d_binding_modes.png")
    compose_2x2(pA, pB, pC, pD, out_fig)
    # Clean temporary single-model PDB conversions
    for f in glob.glob(os.path.join(tmp_dir, "*.pdb")):
        os.remove(f)
    try:
        os.rmdir(tmp_dir)
    except OSError:
        pass
    cmd.quit()


if __name__ == "__main__":
    main()
