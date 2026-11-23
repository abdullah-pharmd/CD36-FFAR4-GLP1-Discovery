"""
02_generate_graphical_abstract.py
Generates a clean, high-resolution, minimal-text Graphical Abstract (300 DPI PNG + LZW TIFF)
formatted for Springer Nature (Molecular Diversity) combining:
  - Left: 2D ACS-clean chemical depictions of flagship Tier 1A leads (Naringenin & Linagliptin/Teneligliptin)
  - Center: Ray-traced PyMOL 3D views of 4A5S (DPP-4 pocket) & 5LGD (CD36 hydrophobic tunnel)
  - Right: Concise quantitative validation badges (N=60 ROC-AUC, 0.78 Å RMSD, PROPKA 3.5, OpenMM MD)
"""
import os
import json
import pandas as pd
from PIL import Image, ImageDraw, ImageFont, ImageChops
from rdkit import Chem
from rdkit.Chem import Draw

SUB_DIR = os.path.dirname(os.path.abspath(__file__))
PROT_B_DIR = os.path.dirname(SUB_DIR)
ROOT_DIR = os.path.dirname(PROT_B_DIR)
RESULTS_DIR = os.path.join(PROT_B_DIR, "results")
NET_DIR = os.path.join(RESULTS_DIR, "network")
PYMOL_DIR = os.path.join(PROT_B_DIR, "figures", "pymol_3d")


def render_mol_img(smiles, size=(420, 260)):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol, size=size, kekulize=True)
    return img.convert("RGB")


def fit_pymol_panel_preserve_aspect(png_path, box_w, box_h, pad=18):
    """Trims outer pure-white margins and scales proportionally onto a clean white box without distortion."""
    im = Image.open(png_path).convert("RGB")
    bg = Image.new("RGB", im.size, (255, 255, 255))
    diff = ImageChops.difference(im, bg)
    bbox = diff.getbbox()
    if bbox:
        x0 = max(0, bbox[0] - pad)
        y0 = max(0, bbox[1] - pad)
        x1 = min(im.width, bbox[2] + pad)
        y1 = min(im.height, bbox[3] + pad)
        im = im.crop((x0, y0, x1, y1))
    scale = min((box_w - 12) / im.width, (box_h - 12) / im.height)
    nw, nh = max(1, int(round(im.width * scale))), max(1, int(round(im.height * scale)))
    im_rs = im.resize((nw, nh), Image.Resampling.LANCZOS)
    frame = Image.new("RGB", (box_w, box_h), (255, 255, 255))
    # Shift slightly down so top-left badge never touches labels
    ox = (box_w - nw) // 2
    oy = max(4, (box_h - nh) // 2 + 8)
    oy = min(oy, box_h - nh)
    frame.paste(im_rs, (ox, oy))
    return frame


def main():
    df = pd.read_csv(os.path.join(RESULTS_DIR, "protocol_b_cd36_dpp4_complete_results.csv")).set_index("Compound")
    df_conv = pd.read_csv(os.path.join(RESULTS_DIR, "protocol_b_convergence_and_rmsd.csv"))
    with open(os.path.join(RESULTS_DIR, "redocking_metrics.json")) as f:
        redock = json.load(f)
    with open(os.path.join(ROOT_DIR, "results", "roc_summary_all.json")) as f:
        roc = json.load(f)["DPP4"]
    with open(os.path.join(NET_DIR, "empirical_target_gap_report.json")) as f:
        gap = json.load(f)
    with open(os.path.join(RESULTS_DIR, "propka_pka_audit_pH74.json")) as f:
        pjson = json.load(f)

    n_roc = int(roc["N_Actives"]) + int(roc["N_Inactives"])
    leads_conv = df_conv[~df_conv["Compound"].str.contains("PLM-512", na=False)]
    max_lead_drift = float(leads_conv["Conv_Pose_RMSD_A"].max())
    cd36_flips = int(pjson["5LGD_CD36"]["Pocket_Focus_Flips_at_pH74"])
    cd36_tot = int(pjson["5LGD_CD36"]["Pocket_Focus_Count"])
    dpp4_flips = int(pjson["4A5S_DPP4"]["Pocket_Focus_Flips_at_pH74"])
    dpp4_tot = int(pjson["4A5S_DPP4"]["Pocket_Focus_Count"])
    emp_nodes = int(gap["Empirical_STRING_Stats"]["total_nodes"])
    dpp4_hits = len(gap["DPP4_Compound_Hits"])
    cd36_hits = len(gap["CD36_Compound_Hits"])

    # Check if Colab MD summary is available
    md_csv = os.path.join(PROT_B_DIR, "colab_md", "results", "md_summary_all_complexes.csv")
    has_md = os.path.exists(md_csv)
    md_ns = int(round(float(pd.read_csv(md_csv)["Production_ns"].iloc[0]))) if has_md else 0

    # Canvas: 2400 x 1050 px @ 300 DPI (8.0 x 3.5 in — ideal TOC graphical abstract ratio)
    W, H = 2400, 1050
    canvas = Image.new("RGB", (W, H), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)

    try:
        f_hdr = ImageFont.truetype("arialbd.ttf", 36)
        f_sub = ImageFont.truetype("arialbd.ttf", 26)
        f_lbl = ImageFont.truetype("arialbd.ttf", 23)
        f_sm = ImageFont.truetype("arial.ttf", 20)
    except Exception:
        f_hdr = f_sub = f_lbl = f_sm = ImageFont.load_default()

    # Top title banner
    draw.rectangle([0, 0, W, 92], fill=(11, 60, 93))
    draw.text((40, 24), "Dual Hepatic CD36 (5LGD) & Vascular DPP-4 (4A5S) Prioritization at Plasma pH 7.4", fill=(255, 255, 255), font=f_hdr)

    # --- Zone 1 (Left): Tier 1A Flagship 2D Scaffolds ---
    z1_x0, z1_x1 = 30, 680
    draw.rounded_rectangle([z1_x0, 118, z1_x1, H - 30], radius=18, fill=(248, 250, 252), outline=(2, 132, 199), width=3)
    draw.text((z1_x0 + 24, 136), "Tier 1A Dual-Target Leads (0 PAINS)", fill=(3, 105, 161), font=f_sub)

    nar_im = render_mol_img(df.loc["Naringenin", "SMILES"], (590, 270))
    ten_im = render_mol_img(df.loc["Teneligliptin", "SMILES"], (590, 270))
    canvas.paste(nar_im, (z1_x0 + 30, 190))
    draw.text((z1_x0 + 42, 472), f"Naringenin (Tier 1A-NP, LogS {df.loc['Naringenin','ESOL_LogS_mol_L']:.2f})", fill=(21, 128, 61), font=f_lbl)
    draw.text((z1_x0 + 42, 502), f"DPP-4 {df.loc['Naringenin','DPP4_4A5S_Affinity_kcal_mol']:.2f} | CD36 {df.loc['Naringenin','CD36_5LGD_Affinity_kcal_mol']:.2f} kcal/mol", fill=(51, 65, 85), font=f_sm)

    draw.line([(z1_x0 + 30, 545), (z1_x1 - 30, 545)], fill=(203, 213, 225), width=2)
    canvas.paste(ten_im, (z1_x0 + 30, 565))
    draw.text((z1_x0 + 42, 845), "Teneligliptin (Tier 1A-Gliptin, Class III)", fill=(15, 118, 110), font=f_lbl)
    draw.text((z1_x0 + 42, 875), f"DPP-4 {df.loc['Teneligliptin','DPP4_4A5S_Affinity_kcal_mol']:.2f} | CD36 {df.loc['Teneligliptin','CD36_5LGD_Affinity_kcal_mol']:.2f} kcal/mol", fill=(51, 65, 85), font=f_sm)
    draw.text(
        (z1_x0 + 42, 914),
        f"Linagliptin ({df.loc['Linagliptin','DPP4_4A5S_Affinity_kcal_mol']:.2f}/{df.loc['Linagliptin','CD36_5LGD_Affinity_kcal_mol']:.2f}) & Sitagliptin ({df.loc['Sitagliptin','DPP4_4A5S_Affinity_kcal_mol']:.2f}/{df.loc['Sitagliptin','CD36_5LGD_Affinity_kcal_mol']:.2f})",
        fill=(71, 85, 105), font=f_sm,
    )

    # Arrow 1 -> 2
    draw.line([(692, 560), (748, 560)], fill=(11, 60, 93), width=8)
    draw.polygon([(748, 542), (776, 560), (748, 578)], fill=(11, 60, 93))

    # --- Zone 2 (Center): PyMOL 3D Ray-Traced Pockets (1:1 Aspect Ratio) ---
    z2_x0, z2_x1 = 788, 1680
    draw.rounded_rectangle([z2_x0, 118, z2_x1, H - 30], radius=18, fill=(255, 255, 255), outline=(15, 118, 110), width=3)
    draw.text((z2_x0 + 24, 136), "3D Ray-Traced Binding Modes (pH 7.4, PROPKA Verified)", fill=(15, 118, 110), font=f_sub)

    p_dpp4 = fit_pymol_panel_preserve_aspect(os.path.join(PYMOL_DIR, "panelB_dpp4_tier1_leads.png"), 846, 375)
    p_cd36 = fit_pymol_panel_preserve_aspect(os.path.join(PYMOL_DIR, "panelD_cd36_tunnel_contacts.png"), 846, 375)
    canvas.paste(p_dpp4, (z2_x0 + 23, 186))
    draw.rectangle([z2_x0 + 23, 186, z2_x0 + 869, 561], outline=(203, 213, 225), width=2)
    draw.rectangle([z2_x0 + 33, 196, z2_x0 + 430, 234], fill=(255, 255, 255), outline=(15, 118, 110), width=2)
    draw.text((z2_x0 + 44, 202), "DPP-4 (4A5S) S1/S2/S2' Pocket", fill=(11, 60, 93), font=f_lbl)

    canvas.paste(p_cd36, (z2_x0 + 23, 590))
    draw.rectangle([z2_x0 + 23, 590, z2_x0 + 869, 965], outline=(203, 213, 225), width=2)
    draw.rectangle([z2_x0 + 33, 600, z2_x0 + 505, 638], fill=(255, 255, 255), outline=(185, 28, 28), width=2)
    draw.text((z2_x0 + 44, 606), "CD36 (5LGD) Segmental Lipid Tunnel", fill=(153, 27, 27), font=f_lbl)

    # Arrow 2 -> 3
    draw.line([(1692, 560), (1748, 560)], fill=(11, 60, 93), width=8)
    draw.polygon([(1748, 542), (1776, 560), (1748, 578)], fill=(11, 60, 93))

    # --- Zone 3 (Right): Falsifiable Benchmarks & Therapeutic Hypothesis ---
    z3_x0, z3_x1 = 1788, W - 30
    draw.rounded_rectangle([z3_x0, 118, z3_x1, H - 30], radius=18, fill=(254, 252, 232), outline=(217, 119, 6), width=3)
    draw.text((z3_x0 + 24, 136), "Validation & Output", fill=(146, 64, 14), font=f_sub)

    cards = [
        ("Empirical Target Gap", f"{emp_nodes} unspiked STRING hubs\nDPP4 = {dpp4_hits} hits | CD36 = {cd36_hits} hits", (240, 249, 255), (2, 132, 199)),
        ("Benchmark Fidelity", f"N7F-901 RMSD = {redock['N7F_901_Hungarian_RMSD_A']:.3f} Å\nN={n_roc} ROC-AUC = {roc['ROC_AUC_Affinity']:.3f} (p={roc['MannWhitney_p_Affinity']:.4f})", (240, 253, 244), (22, 163, 74)),
        ("Electrostatics & Physics", f"PROPKA 3.5: {cd36_flips}/{cd36_tot} CD36, {dpp4_flips}/{dpp4_tot} DPP-4\n" + (f"{md_ns} ns OpenMM MD + MM-GBSA" if has_md else f"Lead exh=8→16 drift ≤ {max_lead_drift:.2f} Å"), (245, 243, 255), (124, 58, 237)),
        ("Dual Hypothesis (MASLD)", "↓ Palmitate flux (CD36)\n↑ Endogenous GLP-1 (DPP-4)\n→ Ready for IC50 & SPR/ITC", (254, 242, 242), (220, 38, 38)),
    ]
    cy = 192
    for title_c, body_c, fc, ec in cards:
        h_c = 186 if "\n\n" in body_c or body_c.count("\n") >= 2 else 166
        draw.rounded_rectangle([z3_x0 + 22, cy, z3_x1 - 22, cy + h_c], radius=12, fill=fc, outline=ec, width=2)
        draw.text((z3_x0 + 38, cy + 16), title_c, fill=ec, font=f_lbl)
        draw.multiline_text((z3_x0 + 38, cy + 52), body_c, fill=(30, 41, 59), font=f_sm, spacing=8)
        cy += h_c + 24

    out_png = os.path.join(SUB_DIR, "Graphical_Abstract_300dpi.png")
    out_tif = os.path.join(SUB_DIR, "Graphical_Abstract_300dpi.tiff")
    canvas.save(out_png, dpi=(300, 300))
    canvas.save(out_tif, dpi=(300, 300), compression="tiff_lzw")
    print(f"[SUCCESS] Saved Graphical Abstract: {out_png} ({os.path.getsize(out_png)/1024:.1f} KB) and {out_tif} ({os.path.getsize(out_tif)/1024:.1f} KB)")


if __name__ == "__main__":
    main()
