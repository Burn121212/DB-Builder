# ==================== ATLASMX ITS — DB BUILDER V9.1 FULL PIPELINE (FOR_R EXPORT) ====================
#
#
# PART 3 (For_R export):
# 12) Create For_R/ folder with phyloseq-ready files for each final DB produced:
#     For_R/<tag>/taxonomy.csv
#     For_R/<tag>/abundance_table.csv
#     For_R/<tag>/sequences.fasta
#     For_R/<tag>/fungal_traits_table.csv
#
# NOTE:
# - MODE controls which tables are used for collapse in Part 2.
# - Regardless of MODE, if MODE="all_eukaryotes" we also export a fungi-only final DB + its For_R package.
# - For_R folder names are strategy-specific and threshold-specific to avoid overwriting.

import pandas as pd
from pathlib import Path
from collections import defaultdict
import time
import re
import sys

# ==================== Check arguments ====================

try:
    OUT_PREFIX_name = sys.argv[1]
    MODE_name = sys.argv[2]
    COLLAPSE_STRATEGY_name = sys.argv[3]
    P_VALUE_THRESHOLD_in = sys.argv[4]
except IndexError:
    print("""
ATLASMX ITS — DB BUILDER V9.1 BUILD R

Usage:
dbbuilder_build_r.py <input_directory> <output_prefix> <collapse_mode> <collapse_strategy> <p-value_threshold>

Output prefix:
    string indicating the path for the DBBuilder output, for example "~/atlas/dbbuilder/db_test"
Collapse mode:
   "fungi"          -> collapse only fungi tables
   "all_eukaryotes" -> collapse full unfiltered tables (all eukaryotes)
Collapse strategy:
   "species_only" -> collapse only when species exists and species_pvalue >= P_VALUE_THRESHOLD
   "genus"        -> collapse only when genus exists and genus_pvalue >= P_VALUE_THRESHOLD
   "all"          -> recursive lowest-rank collapse:
                     species -> genus -> family -> order -> class -> phylum -> domain
                     using the deepest available rank with p >= P_VALUE_THRESHOLD
p-value threshold:
    1.00    ->    very strict
    0.00    ->    similar to vsearch OTU 97%
    
Example:
python snakes/dbbuilder_build_r.py ~/atlas/dbbuilder/db_test all_eukaryotes species_only 1.0
    """)
    sys.exit()

# ====================
# CONFIGURATION
# ====================         

# Collapse mode:
#   "fungi"          -> collapse only fungi tables
#   "all_eukaryotes" -> collapse full unfiltered tables (all eukaryotes)
MODE = MODE_name

# Part 1 I/O
OUT_PREFIX = Path(OUT_PREFIX_name) 
OUTPUT_DIR = OUT_PREFIX / "concatenated_tables"
FUNGI_DIR = OUTPUT_DIR / "Fungi_concatenated"
NON_ANNOTATED_DIR = OUTPUT_DIR / "Non_annotated"
ALL_EUK_DIR = OUTPUT_DIR / "AllEuk_concatenated"

# Part 2 I/O
FINAL_DB_DIR = OUT_PREFIX / "FINAL_DB" 
FUNGAL_TRAITS_DB = Path("snakes/Fungal_Traits_DB.txt") 

# Part 3 I/O
FOR_R_DIR = OUT_PREFIX / "For_R"

# Create required directories
for directory in [FOR_R_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

# Collapse strategy:
#   "species_only" -> collapse only when species exists and species_pvalue >= P_VALUE_THRESHOLD
#   "genus"        -> collapse only when genus exists and genus_pvalue >= P_VALUE_THRESHOLD
#   "all"          -> recursive lowest-rank collapse:
#                     species -> genus -> family -> order -> class -> phylum -> domain
#                     using the deepest available rank with p >= P_VALUE_THRESHOLD
COLLAPSE_STRATEGY = COLLAPSE_STRATEGY_name

# Main collapse confidence threshold
P_VALUE_THRESHOLD = float(P_VALUE_THRESHOLD_in)    # p= 1.00 very strict  p= 0.00 similar to vsearch OTU 97%

# Threshold used for:
# 1) selecting the deepest valid taxon for SPPN assignment
# 2) sanitizing parsed taxonomy columns (blank low-confidence names)
SPPN_P_THRESHOLD = 0.80


# =========================================================
# ==================== HELPER NAMING FUNCTIONS ============
# =========================================================

def format_threshold_for_name(value: float) -> str:
    """Convert float threshold to filename-safe format, e.g. 1.0 -> p1p0."""
    return f"p{value}".replace(".", "p")


def build_final_db_filename(suffix: str) -> str:
    """
    Build final DB filename using:
    strategy + suffix + collapse threshold
    """
    return f"Final_Database_{COLLAPSE_STRATEGY}_{suffix}_p{P_VALUE_THRESHOLD}.csv"


def build_for_r_tag(suffix: str) -> str:
    """
    Build unique For_R folder name using:
    dataset suffix + strategy + collapse threshold + sppn threshold

    Examples:
      fungi_species_only_p1p0_sppn0p8
      fungi_genus_p1p0_sppn0p8
      fungi_all_p1p0_sppn0p8
      all_eukaryotes_all_p1p0_sppn0p8
    """
    collapse_tag = format_threshold_for_name(P_VALUE_THRESHOLD)
    sppn_tag = f"sppn{SPPN_P_THRESHOLD}".replace(".", "p")
    return f"{suffix}_{COLLAPSE_STRATEGY}_{collapse_tag}_{sppn_tag}"




def load_fungal_traits_db():
    """Load Fungal Traits DB and normalize its key columns."""
    print("\n📁 LOADING FUNGAL TRAITS DATABASE...")
    if not FUNGAL_TRAITS_DB.exists():
        print(f"   ⚠️ Warning: {FUNGAL_TRAITS_DB} not found. Lifestyle columns will be empty.")
        return None

    try:
        traits = pd.read_csv(FUNGAL_TRAITS_DB, sep="\t", low_memory=False)
        cols = traits.columns.tolist()

        genus_col = None
        primary_col = None
        secondary_col = None

        for c in cols:
            cl = c.lower()
            if cl == "genus":
                genus_col = c
            elif cl == "primary_lifestyle":
                primary_col = c
            elif cl == "secondary_lifestyle":
                secondary_col = c

        if not genus_col:
            genus_col = next((c for c in cols if "genus" in c.lower()), None)
        if not primary_col:
            primary_col = next((c for c in cols if "primary" in c.lower() and "lifestyle" in c.lower()), None)
        if not secondary_col:
            secondary_col = next((c for c in cols if "secondary" in c.lower() and "lifestyle" in c.lower()), None)

        if not all([genus_col, primary_col, secondary_col]):
            print("   ⚠️ Warning: Missing required columns in Fungal Traits DB.")
            return None

        traits = traits[[genus_col, primary_col, secondary_col]].copy()
        traits = traits.rename(columns={
            genus_col: "GENUS",
            primary_col: "primary_lifestyle",
            secondary_col: "Secondary_lifestyle"
        })

        traits["GENUS"] = traits["GENUS"].astype(str)
        traits["GENUS_clean"] = traits["GENUS"].str.strip()
        traits = traits[traits["GENUS_clean"] != ""].copy()

        before = len(traits)
        traits = traits.drop_duplicates(subset=["GENUS_clean"], keep="first").copy()
        after = len(traits)
        if before != after:
            print(f"   ⚠️ Deduplicated traits by GENUS_clean: {before} -> {after}")

        traits["primary_lifestyle"] = traits["primary_lifestyle"].fillna("")
        traits["Secondary_lifestyle"] = traits["Secondary_lifestyle"].fillna("")
        print(f"   ✅ Fungal Traits records (unique GENUS_clean): {len(traits)}")
        return traits

    except Exception as e:
        print(f"   ❌ Error loading fungal traits database: {e}")
        return None
        
# =========================================================
# ==================== PART 3 (FOR_R) ======================
# =========================================================

TAX_RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
PVALUE_COLS = [f"{r}_pvalue" for r in TAX_RANKS]


def detect_site_columns(df):
    """Detect abundance/site columns for phyloseq export."""
    non_site = set([
        "SPPN", "OTU_XX", "Original_OTUID",
        "sequence", "sequence_lenght", "sequence_length",
        "Total_Abundance",
        "OTU_Number", "sintax_taxonomy",
        "primary_lifestyle", "Secondary_lifestyle"
    ] + TAX_RANKS + PVALUE_COLS)

    extra_non_site = [c for c in df.columns if c.lower().startswith("sequence")]
    non_site.update(extra_non_site)

    site_cols = [c for c in df.columns if c not in non_site]

    if not site_cols:
        return site_cols, pd.DataFrame()

    site_df = df[site_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    if ((site_df % 1) != 0).any().any():
        print("⚠️ Warning: some site columns have non-integer values. They will be truncated to integers.")
    site_df = site_df.astype("Int64")
    return site_cols, site_df


def build_taxonomy_csv_for_r(df, outdir: Path):
    """Export taxonomy.csv with SPPN and taxonomy columns."""
    print(f"🔬 Building taxonomy.csv in {outdir} ...")
    missing = [r for r in TAX_RANKS if r not in df.columns]
    if missing:
        raise ValueError(f"Missing taxonomy columns: {missing}")

    tax = df[["SPPN"] + TAX_RANKS].copy().sort_values("SPPN")
    tax.to_csv(outdir / "taxonomy.csv", index=False)
    print(f"   ✅ taxonomy.csv -> {len(tax)} rows")
    return tax


def build_abundance_table_csv_for_r(df, outdir: Path):
    """Export abundance_table.csv with SPPN and abundance counts."""
    print(f"📊 Building abundance_table.csv in {outdir} ...")
    site_cols, site_df = detect_site_columns(df)

    if not site_cols:
        ab = df[["SPPN"]].copy()
        print("   ⚠️ No site columns detected; exporting only SPPN.")
    else:
        ab = pd.concat([df[["SPPN"]].reset_index(drop=True), site_df.reset_index(drop=True)], axis=1)

    ab = ab.sort_values("SPPN")
    ab.to_csv(outdir / "abundance_table.csv", index=False)
    print(f"   ✅ abundance_table.csv -> {len(ab)} rows, {len(ab.columns)-1} sites")
    return ab


def build_sequences_fasta_for_r(df, outdir: Path):
    """Export sequences.fasta using SPPN as FASTA header."""
    print(f"🧬 Building sequences.fasta in {outdir} ...")
    if "sequence" not in df.columns:
        (outdir / "sequences.fasta").write_text("", encoding="utf-8")
        print("   ⚠️ No sequence column; wrote empty FASTA.")
        return 0

    seq_df = df[["SPPN", "sequence"]].dropna(subset=["sequence"]).copy()
    seq_df["sequence"] = seq_df["sequence"].astype(str).str.replace(r"\s+", "", regex=True)
    seq_df = seq_df[seq_df["sequence"].str.len() > 0]

    lines = []
    for _, r in seq_df.iterrows():
        lines.append(f">{r['SPPN']}")
        lines.append(r["sequence"])

    with open(outdir / "sequences.fasta", "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print(f"   ✅ sequences.fasta -> {len(seq_df)} sequences")
    return len(seq_df)


def build_fungal_traits_table_for_r(taxonomy_df, traits_df, outdir: Path):
    """Export fungal_traits_table.csv by joining exact trimmed genus names."""
    print(f"🔍 Building fungal_traits_table.csv in {outdir} ...")
    base = taxonomy_df[["SPPN", "genus"]].copy()
    base["GENUS_clean"] = base["genus"].astype(str).str.strip()

    if traits_df is None or traits_df.empty or "GENUS_clean" not in traits_df.columns:
        out = outdir / "fungal_traits_table.csv"
        base[["SPPN"]].to_csv(out, index=False)
        print(f"   ⚠️ No valid traits DB; exported only SPPN ({len(base)} rows)")
        return base[["SPPN"]]

    merged = base.merge(traits_df, on="GENUS_clean", how="left", suffixes=("", "_traits"))

    matches = merged["GENUS_clean"].isin(traits_df["GENUS_clean"]).sum()
    total = len(merged)
    print(f"   ✅ Exact genus matches: {matches}/{total} ({(matches/total*100 if total else 0):.1f}%)")

    preferred = ["SPPN"]
    for col in ["COMMENT on genus", "primary_lifestyle", "Secondary_lifestyle"]:
        if col in merged.columns:
            preferred.append(col)

    technical = ["genus", "GENUS", "GENUS_clean"]
    other_cols = [c for c in merged.columns if c not in preferred + technical and c != "SPPN"]
    final_cols = preferred + other_cols

    merged = merged[final_cols].drop_duplicates(subset=["SPPN"], keep="first")
    merged.to_csv(outdir / "fungal_traits_table.csv", index=False)
    print(f"   ✅ fungal_traits_table.csv -> {len(merged)} rows")
    return merged


def export_for_r_package(final_db_csv: Path, suffix: str, traits_df):
    """
    Build a phyloseq-ready package under For_R/<tag>/ from a final DB CSV.
    The tag is generated centrally to avoid overwriting between strategies.
    """
    tag = build_for_r_tag(suffix)
    outdir = FOR_R_DIR / tag
    outdir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 70)
    print(f"📦 FOR_R EXPORT: {tag}")
    print("=" * 70)

    df = pd.read_csv(final_db_csv)

    if "SPPN" not in df.columns:
        raise ValueError(f"{final_db_csv} does not contain 'SPPN' column")

    if df["SPPN"].duplicated().any():
        dups = df[df["SPPN"].duplicated(keep=False)]["SPPN"].unique()
        print(f"⚠️ Duplicated SPPN found ({len(dups)}). Keeping first occurrence.")
        df = df.drop_duplicates(subset=["SPPN"], keep="first")

    taxonomy_df = build_taxonomy_csv_for_r(df, outdir)
    abundance_df = build_abundance_table_csv_for_r(df, outdir)
    nseq = build_sequences_fasta_for_r(df, outdir)
    fungal_traits_df = build_fungal_traits_table_for_r(taxonomy_df, traits_df, outdir)

    print("\n📦 SUMMARY (For_R):")
    print(f"   • taxonomy.csv: {len(taxonomy_df)} rows")
    print(f"   • abundance_table.csv: {len(abundance_df)} rows, {len(abundance_df.columns)-1} sites")
    print(f"   • sequences.fasta: {nseq} sequences")
    print(f"   • fungal_traits_table.csv: {len(fungal_traits_df)} rows")
    print(f"   ✅ Ready in: {outdir.resolve()}")


# =========================================================
# ==================== EXPORT FINAL OUTPUTS =================
# =========================================================

def save_outputs(final_df, logs, site_cols, part2_input_dir, part2_pattern):
    """
    Save final DB CSV(s) and the collapse log.
    Filenames include strategy and p-value to avoid overwriting.
    Returns a dict with produced CSV paths for downstream For_R export.
    """
    produced = {}

    suffix_main = "all_eukaryotes" if MODE == "all_eukaryotes" else "fungi"
    out_csv_main = FINAL_DB_DIR / build_final_db_filename(suffix_main)

    rank_cols = TAX_RANKS
    pval_cols = [f"{r}_pvalue" for r in rank_cols]
    lifestyle_cols = ["primary_lifestyle", "Secondary_lifestyle"]

    first_cols = (
        ["SPPN", "OTU_XX", "Original_OTUID"]
        + rank_cols
        + pval_cols
        + lifestyle_cols
        + ["sequence", "sequence_lenght", "Total_Abundance"]
    )
    existing = [c for c in first_cols if c in final_df.columns]
    final_order = existing + [c for c in site_cols if c in final_df.columns]

    final_df[final_order].to_csv(out_csv_main, index=False, encoding="utf-8")
    print(f"💾 Final database saved: {out_csv_main}")
    produced[suffix_main] = out_csv_main

    if MODE == "all_eukaryotes":
        fungi_out = FINAL_DB_DIR / build_final_db_filename("fungi")
        fungi_mask = final_df["domain"].astype(str).str.strip().str.lower().eq("fungi")
        fungi_df = final_df.loc[fungi_mask].copy()
        fungi_df[final_order].to_csv(fungi_out, index=False, encoding="utf-8")
        print(f"🍄 Fungi-only subset saved: {fungi_out} | Rows: {len(fungi_df)}")
        produced["fungi"] = fungi_out

    log_file = FINAL_DB_DIR / build_collapse_log_filename(suffix_main)
    with open(log_file, "w", encoding="utf-8") as f:
        f.write(f"LOG — {COLLAPSE_STRATEGY} collapse (p≥{P_VALUE_THRESHOLD})\n")
        f.write(f"MODE: {MODE}\n")
        f.write(f"COLLAPSE_STRATEGY: {COLLAPSE_STRATEGY}\n")
        f.write(f"PART2_INPUT_DIR: {part2_input_dir}\n")
        f.write(f"PART2_PATTERN: {part2_pattern}\n")
        f.write("=" * 70 + "\n\n")
        f.write("Representative chosen ONLY among rows with COMPLETE taxonomy (d,p,c,o,f,g,s).\n")
        f.write("Taxonomy strings cleaned to remove artifacts like '(Fungi)'.\n")
        f.write(f"Parsed taxonomy columns sanitized: names with p<{SPPN_P_THRESHOLD} blanked (cascade).\n")
        f.write("Lifestyle info merged from Fungal_Traits_DB by genus (after sanitization).\n")
        f.write("Rows with Total_Abundance==0 were pruned.\n")
        if MODE == "all_eukaryotes":
            f.write("EXTRA: Exported fungi-only subset from all-eukaryotes final DB.\n")
        f.write("\n")

        removed_total = 0
        collapsed_groups = 0
        skipped_groups = 0

        success_note = f"collapsed_with_complete_tax_{COLLAPSE_STRATEGY}"
        skip_note = f"skip_collapse_incomplete_tax_{COLLAPSE_STRATEGY}"

        for i, L in enumerate(logs, 1):
            removed_total += L.get("removed_count", 0)
            if L.get("note") == success_note:
                collapsed_groups += 1
            elif L.get("note") == skip_note:
                skipped_groups += 1

            f.write(f"[Group {i}] key: {L.get('collapse_key')}\n")
            f.write(f"  Representative: {L.get('kept_otu')}\n")
            f.write(f"  Representative length: {L.get('representative_seq_len')}\n")
            f.write(f"  Representative taxonomy: {L.get('representative_tax')}\n")
            f.write(f"  Members: {L.get('group_size')} | Removed: {L.get('removed_count')}\n")
            if L.get("removed_otus"):
                f.write(f"  Removed (up to 10): {', '.join(L['removed_otus'])}\n")
            f.write(f"  Note: {L.get('note')}\n\n")

        f.write("SUMMARY\n")
        f.write(f"  Total removed: {removed_total}\n")
        f.write(f"  Final OTUs (after pruning zeros): {len(final_df)}\n")
        f.write(f"  Groups successfully collapsed: {collapsed_groups}\n")
        f.write(f"  Groups skipped (incomplete taxonomy): {skipped_groups}\n")

    print(f"📋 Collapse log saved: {log_file}")
    return produced


# =========================================================
# ==================== MAIN (PART 1 + PART 2 + PART 3) =====
# =========================================================

def main():
    print("=" * 70)
    print("🚀 ATLASMX ITS — DB BUILDER V9.1")
    print("=" * 70)
    print("PART 3: For_R phyloseq-ready export packages")
    print(f"MODE: {MODE} | strategy: {COLLAPSE_STRATEGY} | collapse p≥{P_VALUE_THRESHOLD} | sanitize p≥{SPPN_P_THRESHOLD}")
    print("=" * 70)
    print()

    # load fungal traits
    fungal_traits = load_fungal_traits_db()
    
    # get file names from collapse step
    produced_csvs = {}
    suffix_main = "all_eukaryotes" if MODE == "all_eukaryotes" else "fungi"
    out_csv_main = FINAL_DB_DIR / build_final_db_filename(suffix_main)
    produced_csvs[suffix_main] = out_csv_main
    if MODE == "all_eukaryotes":
        fungi_out = FINAL_DB_DIR / build_final_db_filename("fungi")
        produced_csvs["fungi"] = fungi_out
    
    # ---------------- PART 3 ----------------
    print("\n🧬 PART 3 — EXPORTING For_R PACKAGES")    

    if "all_eukaryotes" in produced_csvs:
        export_for_r_package(
            produced_csvs["all_eukaryotes"],
            suffix="all_eukaryotes",
            traits_df=fungal_traits
        )

    if "fungi" in produced_csvs:
        export_for_r_package(
            produced_csvs["fungi"],
            suffix="fungi",
            traits_df=fungal_traits
        )

    print("\n🎉 ALL DONE")
    print(f"   MODE: {MODE}")
    print(f"   COLLAPSE_STRATEGY: {COLLAPSE_STRATEGY}")
    print(f"   For_R folder: {FOR_R_DIR.resolve()}")


if __name__ == "__main__":
    main()