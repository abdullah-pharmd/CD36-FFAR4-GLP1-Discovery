"""
01_stylometric_and_ai_tell_audit.py
Pre-submission editorial, stylometric, and AI-pattern hygiene auditor for Springer Nature
(Molecular Diversity) manuscripts. Checks:
  1. Banned LLM cliché / buzzword density (strict zero-tolerance list + connector caps)
  2. Sentence-length burstiness (mean, SD, CV >= 0.35, short/long sentence balance)
  3. Em-dash overuse, repeated paragraph openers, and non-technical 4-gram repetition
  4. Springer Nature structural compliance (Abstract <= 250 words & 0 citations, 4-6 Keywords,
     mandatory Declarations subheadings, sequential [1]..[N] citations, 0 raw formatting leaks)
"""
import os
import re
import sys
import math
from collections import Counter
import docx

BANNED_BUZZWORDS = [
    "delve", "delves", "delving", "crucially", "critically,", "underscore", "underscores",
    "underscoring", "pivotal", "tapestry", "realm", "multifaceted", "holistic", "paramount",
    "meticulously", "orchestrate", "orchestrates", "orchestrating", "paradigm shift",
    "sheds light", "shed light", "paves the way", "pave the way", "synthesizing",
    "leveraging", "harnessing", "unprecedented", "groundbreaking", "game-changer",
    "cornerstone", "beacon", "intricate interplay", "synergistic interplay", "testament",
    "bustling", "vibrant", "ever-evolving", "fostering", "cutting-edge", "state-of-the-art",
]

CAPPED_CONNECTORS = {
    "notably": 2,
    "importantly": 2,
    "furthermore": 2,
    "moreover": 2,
    "additionally": 2,
    "consequently": 3,
    "specifically": 4,
    "however": 8,
}

REQUIRED_DECLARATIONS = [
    "Funding",
    "Competing Interests",
    "Ethics Approval",
    "Consent to Participate",
    "Consent for Publication",
    "Data and Code Availability",
    "Author Contributions",
]


def split_sentences(text):
    # Protect common scientific abbreviations and decimal numbers
    t = re.sub(r"\b(et al|Fig|Figs|Eq|Table|Dr|vs|i\.e|e\.g|approx|min|max|No|Ref|Refs)\.", r"\1<DOT>", text)
    t = re.sub(r"(\d)\.(\d)", r"\1<DOT>\2", t)
    sents = [s.replace("<DOT>", ".").strip() for s in re.split(r"(?<=[.!?])\s+", t) if len(s.strip()) > 12]
    return sents


def audit_manuscript(docx_path):
    doc = docx.Document(docx_path)
    paras = [p.text.strip() for p in doc.paragraphs if p.text.strip()]
    tables_text = "\n".join(c.text for t in doc.tables for r in t.rows for c in r.cells)

    # Split body before References (and isolate Abstract)
    ref_idx = next((i for i, p in enumerate(paras) if re.match(r"^(References|\d+\.\s*References)$", p, re.I)), len(paras))
    abs_idx = next((i for i, p in enumerate(paras) if p.lower() == "abstract"), -1)
    kw_idx = next((i for i, p in enumerate(paras) if p.lower().startswith("keywords")), -1)
    intro_idx = next((i for i, p in enumerate(paras) if re.match(r"^1\.\s*Introduction", p, re.I)), -1)

    abstract_txt = " ".join(paras[abs_idx + 1:kw_idx]) if (abs_idx >= 0 and kw_idx > abs_idx) else ""
    keywords_txt = paras[kw_idx] if kw_idx >= 0 else ""
    prose_paras = [
        p for p in paras[max(0, intro_idx):ref_idx]
        if not re.match(r"^(Figure|Table)\s+\d+", p)
        and not re.match(r"^\d+(\.\d+)?\s+[A-Z]", p)
        and not p.startswith("[Eq.")
        and "RMSD_Hungarian" not in p
        and "ΔG_bind =" not in p
        and "LE = " not in p
        and len(p.split()) >= 18
    ]
    prose_all = "\n".join(prose_paras)
    body_all = "\n".join(paras[:ref_idx])

    issues = []
    warnings = []

    print("=" * 78)
    print(f"STYLOMETRIC, AI-TELL & SPRINGER COMPLIANCE AUDIT: {os.path.basename(docx_path)}")
    print("=" * 78)

    # 1. Word count & Abstract checks
    total_words = len(prose_all.split()) + len(abstract_txt.split())
    abs_words = len(abstract_txt.split())
    print(f"Total Prose Word Count (Abstract + Body): {total_words:,} words")
    print(f"Abstract Word Count: {abs_words} words (Springer Molecular Diversity target <= 250)")
    if abs_words == 0:
        issues.append("Abstract not found.")
    elif abs_words > 250:
        issues.append(f"Abstract too long ({abs_words} words; Springer Molecular Diversity strict limit <= 250).")
    if re.search(r"\[\d+\]", abstract_txt):
        issues.append("Citation bracket found inside Abstract (Springer prohibits citations in Abstract).")

    kw_clean = re.sub(r"^keywords\s*:\s*", "", keywords_txt, flags=re.I)
    kw_list = [k.strip() for k in re.split(r"[;,]", kw_clean) if k.strip()]
    print(f"Keywords ({len(kw_list)}): {kw_list}")
    if not (4 <= len(kw_list) <= 6):
        issues.append(f"Expected 4-6 keywords, got {len(kw_list)}.")

    # 2. Banned AI Buzzword Scan
    low_prose = (abstract_txt + "\n" + prose_all).lower()
    found_banned = {}
    for bw in BANNED_BUZZWORDS:
        cnt = len(re.findall(r"\b" + re.escape(bw.rstrip(",")) + (r"," if bw.endswith(",") else r"\b"), low_prose))
        if cnt > 0:
            found_banned[bw] = cnt
    print(f"Banned AI Buzzwords Found: {len(found_banned)} {found_banned if found_banned else '[CLEAN - 0]'}")
    if found_banned:
        issues.append(f"Banned AI buzzwords present: {found_banned}")

    # 3. Capped Transition Connector Check
    over_conn = {}
    for conn, cap in CAPPED_CONNECTORS.items():
        cnt = len(re.findall(r"\b" + re.escape(conn) + r"\b", low_prose))
        if cnt > cap:
            over_conn[conn] = f"{cnt} (cap {cap})"
    print(f"Overused Transition Connectors: {over_conn if over_conn else '[CLEAN]'}")
    if over_conn:
        warnings.append(f"Overused transitions: {over_conn}")

    # 4. Sentence-Length Burstiness & Em-Dash Density
    sents = split_sentences(abstract_txt + " " + prose_all)
    lens = [len(s.split()) for s in sents]
    mean_l = sum(lens) / max(1, len(lens))
    sd_l = math.sqrt(sum((x - mean_l) ** 2 for x in lens) / max(1, len(lens)))
    cv = sd_l / mean_l if mean_l > 0 else 0.0
    pct_short = 100.0 * sum(1 for x in lens if x <= 14) / max(1, len(lens))
    pct_long = 100.0 * sum(1 for x in lens if x >= 30) / max(1, len(lens))
    em_dashes = prose_all.count("—")
    em_per_1k = 1000.0 * em_dashes / max(1, total_words)
    print(f"Sentence Rhythm ({len(sents)} sentences): Mean={mean_l:.1f}w | SD={sd_l:.1f}w | Burstiness CV={cv:.3f} (target >= 0.35)")
    print(f"  Short (<=14w): {pct_short:.1f}% | Long (>=30w): {pct_long:.1f}% | Em-dashes: {em_dashes} ({em_per_1k:.1f}/1k words, target <= 8)")
    if cv < 0.35:
        warnings.append(f"Low sentence burstiness CV={cv:.3f} (<0.35); vary sentence lengths more.")
    if mean_l > 27.5:
        warnings.append(f"High mean sentence length ({mean_l:.1f}w > 27.5w); split dense multi-clause sentences.")
    if pct_short < 12.0:
        warnings.append(f"Low short-sentence share ({pct_short:.1f}% < 12%); mix in crisp short sentences.")
    if em_per_1k > 8.0:
        warnings.append(f"High em-dash density ({em_per_1k:.1f}/1k words); replace excess em-dashes with commas or parentheses.")

    # 5. Repeated Paragraph Openers & Internal Repo Jargon Leak Check
    openers = Counter(" ".join(p.split()[:2]).lower() for p in prose_paras)
    rep_openers = {k: v for k, v in openers.items() if v >= 3}
    print(f"Repeated Paragraph Openers (>=3x): {rep_openers if rep_openers else '[CLEAN]'}")
    if rep_openers:
        warnings.append(f"Repeated paragraph openers: {rep_openers}")
    for jargon in ["Protocol B", "Protocol A", "Update 1", "Update 2", "Update 3"]:
        if jargon.lower() in body_all.lower():
            issues.append(f"Internal repository codename '{jargon}' leaked into standalone journal manuscript.")

    # 6. Mandatory Springer Declarations
    decl_missing = [d for d in REQUIRED_DECLARATIONS if d.lower() not in body_all.lower()]
    print(f"Missing Mandatory Springer Declarations: {decl_missing if decl_missing else '[ALL PRESENT]'}")
    if decl_missing:
        issues.append(f"Missing Springer Declarations: {decl_missing}")

    # 7. Formatting Leak & Citation Bijection Check
    full_comb = "\n".join(paras) + "\n" + tables_text
    for bad in ["$", "\\", "**", "(False)", "(True)", "(false)", "(true)", "(0P)", "((", "+/-", "-C#N", "Khan¹*¹"]:
        cnt = full_comb.count(bad)
        if cnt > 0:
            issues.append(f"Forbidden token '{bad}' appeared {cnt}x.")

    # Citation sequential check
    ref_lines = [p for p in paras[ref_idx + 1:] if re.match(r"^\d+\.\s+", p)]
    n_refs = len(ref_lines)
    cited = []
    for grp in re.findall(r"\[(\d+(?:[,\-–]\s*\d+)*)\]", body_all):
        for part in grp.split(","):
            part = part.strip()
            if "-" in part or "–" in part:
                a, b = re.split(r"[-–]", part)
                cited.extend(range(int(a.strip()), int(b.strip()) + 1))
            elif part.isdigit():
                cited.append(int(part))
    seen, seq = set(), []
    for n in cited:
        if n not in seen:
            seen.add(n)
            seq.append(n)
    missing_refs = sorted(set(range(1, n_refs + 1)) - seen)
    extra_refs = sorted(seen - set(range(1, n_refs + 1)))
    print(f"References: {n_refs} entries | Cited unique: {len(seen)} | Sequential: {seq == sorted(seq)} | Missing: {missing_refs} | Extra: {extra_refs}")
    if missing_refs or extra_refs:
        issues.append(f"Citation mismatch: missing={missing_refs}, extra={extra_refs}")
    if seq != sorted(seq):
        issues.append(f"Citations out of sequential order: first-seen={seq}")

    print("-" * 78)
    if issues:
        print(f"[FAIL] {len(issues)} blocking issue(s):")
        for x in issues:
            print("  ERROR:", x)
    else:
        print("[PASS] Zero blocking issues! Stylometry, Springer layout, and integrity checks passed.")
    if warnings:
        print(f"[NOTES] {len(warnings)} style suggestion(s):")
        for w in warnings:
            print("  WARN:", w)
    return len(issues) == 0


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "02_Manuscript_Molecular_Diversity.docx"
    )
    ok = audit_manuscript(target)
    sys.exit(0 if ok else 1)
