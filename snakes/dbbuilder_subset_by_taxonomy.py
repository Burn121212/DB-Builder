# ============================================================
#  SUBSET BY TAXONOMY + PRUNE (OTUs & SITES) + DROP INPUT LIFESTYLES + MERGE FULL FUNGAL TRAITS
#
#  INPUT columns expected (example):
#  SPPN, OTU_XX, Original_OTUID, domain..species, *_pvalue,
#  primary_lifestyle, Secondary_lifestyle, sequence, sequence_lenght, Total_Abundance,
#  Sitio_AtlasEXP01, Sitio_AtlasEXP02, ...
#
#  OUTPUT columns desired:
#  [fixed columns WITHOUT input lifestyles] + [renamed sites] + [FULL fungal traits columns]
#  where primary_lifestyle & Secondary_lifestyle come ONLY from traits.
# ============================================================

from pathlib import Path
import pandas as pd
import re
import sys

# ==================== Check arguments ====================

try:
    OUT_PREFIX_name = sys.argv[1]
    INPUT_File_name = sys.argv[2]
    TAXONOMIC_RANK_name = sys.argv[3]
    TAXONOMIC_FILTER_name = sys.argv[4]
    OUTPUT_name = sys.argv[5]
    SUFFIX_TOTAL_ABUNDANCE_name = sys.argv[6]
    REMOVE_PREFIX_name = sys.argv[7]
    TAX_FILTER_MODE_name = sys.argv[8]
    
except IndexError:
    print("""
ATLASMX ITS — DB BUILDER V9.1
SUBSET BY TAXONOMY + PRUNE (OTUs & SITES) + DROP INPUT LIFESTYLES + MERGE FULL FUNGAL TRAITS

Usage:
dbbuilder_subset_by_taxonomy.py <output_prefix> <input_file_name> <taxonomic_rank> <taxonomic_filter> <output_name> <suffix_total_abundance> <remove_prefix> <tax_filter_mode>

Output prefix:
    string indicating the path of the DBBuilder output, for example "~/atlas/dbbuilder/db_test"
Input file name:
    string for the input file, for example "FINAL_DB/Final_Database_species_only_all_eukaryotes_p0.8.csv"
Taxonomic rank:
    "domain" | "phylum" | "class" | "order" | "family" | "genus" | "species"
Taxonomic filter:
    string, for example "Amanitaceae"
    list, for example "Amanitaceae,Boletaceae"
Output name:
    string for output  file names, for example "MXABC_ITS_unoise_zOTUs_subset_family_Amanitaceae_p1.0_species"
Suffix total abundance:
    string to reame total abundance, for example "MXABC_zOTUs_Amanitaceae_p1.0_species"
Remove prefix:
    string to remove from sample names, for example "Sitio_"
Taxonomic filter mode:
    "exact" | "loose" | "contains"
    
Example:
python snakes/dbbuilder_subset_by_taxonomy.py ~/atlas/dbbuilder/db_test FINAL_DB/Final_Database_species_only_all_eukaryotes_p0.8.csv family Amanitaceae MXABC_ITS_unoise_zOTUs_subset_family_Amanitaceae_p1.0_species MXABC_zOTUs_Amanitaceae_p1.0_species Sitio_ exact
    """)
    sys.exit()


# =====================
# CONFIGURATION 
# ===================== 

OUT_PREFIX = Path(OUT_PREFIX_name) 
INPUT_FILE = OUT_PREFIX / INPUT_File_name #  needs fine tuning 
FUNGAL_TRAITS_DB = Path("snakes/Fungal_Traits_DB.txt")  # TSV
OUTPUT_DIR = OUT_PREFIX / "subset"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Taxonomic filter
TAXONOMIC_RANK = TAXONOMIC_RANK_name

TAXONOMIC_FILTER = TAXONOMIC_FILTER_name  # str or list
if TAXONOMIC_FILTER == 'None':
    TAXONOMIC_FILTER = None
elif ',' in TAXONOMIC_FILTER:   # get list if a comma is present
    TAXONOMIC_FILTER = TAXONOMIC_FILTER.split(',')

TAX_FILTER_MODE = TAX_FILTER_MODE_name

# Output name
OUTPUT_NAME = OUTPUT_name

# Rename Total_Abundance
SUFFIX_TOTAL_ABUNDANCE = SUFFIX_TOTAL_ABUNDANCE_name

# Prune OTUs (rows) with zero reads across ALL sites
DROP_ZERO_OTUS = True

# Prune SITES (columns) whose total sum across all OTUs is zero
DROP_ZERO_SITES = True

# Rename site columns: remove "Sitio_" prefix if present
REMOVE_SITE_PREFIX = REMOVE_PREFIX_name   # set to None to disable
if REMOVE_SITE_PREFIX == 'None':
    REMOVE_SITE_PREFIX = None

# ===================== CONSTANTS =====================

TAX_RANKS = ["domain","phylum","class","order","family","genus","species"]
PVALUE_COLS = [f"{r}_pvalue" for r in TAX_RANKS]

# These are the base columns we always expect / want at the front (order matters)
BASE_FIXED_COLS = [
    "SPPN","OTU_XX","Original_OTUID",
    "domain","phylum","class","order","family","genus","species",
    "domain_pvalue","phylum_pvalue","class_pvalue","order_pvalue","family_pvalue","genus_pvalue","species_pvalue",
    "sequence","sequence_lenght","Total_Abundance"
]

# Columns to remove from input BEFORE merging traits (to avoid duplicates)
DROP_INPUT_LIFESTYLES = {"primary_lifestyle", "Secondary_lifestyle"}

# ===================== HELPERS =====================

def normalize_name(x: str) -> str:
    return str(x).strip().lower()

def apply_taxonomic_filter(df: pd.DataFrame, rank: str, filt, mode: str) -> pd.DataFrame:
    if rank not in df.columns:
        raise ValueError(f"Taxonomic rank '{rank}' not found in columns.")

    if isinstance(filt, (str, int, float)):
        filters = [str(filt)]
    else:
        filters = [str(x) for x in filt]

    col = df[rank].astype(str)

    if mode == "exact":
        mask = col.isin(filters)

    elif mode == "loose":
        filters_norm = {normalize_name(x) for x in filters}
        mask = col.apply(normalize_name).isin(filters_norm)

    elif mode == "contains":
        pattern = "|".join(re.escape(x) for x in filters)
        mask = col.str.contains(pattern, case=False, na=False, regex=True)

    else:
        raise ValueError("TAX_FILTER_MODE must be 'exact', 'loose', or 'contains'.")

    out = df[mask].copy()
    print(f"🔬 Filter applied: {rank} | mode={mode} | values={filters}")
    print(f"   • Rows before: {len(df)}")
    print(f"   • Rows after : {len(out)}")
    return out

def detect_site_columns(df: pd.DataFrame) -> list:
    """
    Detect site columns as everything NOT in fixed/tax/pvalue/known non-site columns.
    We purposely include input lifestyles as non-site so we NEVER parse them as numeric.
    """
    non_site = set(BASE_FIXED_COLS) | set(TAX_RANKS) | set(PVALUE_COLS) | DROP_INPUT_LIFESTYLES | {
        "sintax_taxonomy", "OTU_Number"
    }

    # also protect any sequence* columns
    for c in df.columns:
        if str(c).strip().lower().startswith("sequence"):
            non_site.add(c)

    sites = [c for c in df.columns if c not in non_site]
    return sites

def rename_sites(cols: list, prefix: str | None) -> list:
    if not prefix:
        return cols
    out = []
    for c in cols:
        if isinstance(c, str) and c.startswith(prefix):
            out.append(c[len(prefix):])
        else:
            out.append(c)
    return out

def load_traits(tsv_path: Path) -> pd.DataFrame:
    """
    Load full traits table (keep all columns).
    Must contain a genus column (GENUS/genus/etc).
    """
    if not tsv_path.exists():
        print(f"⚠️ Traits file not found: {tsv_path} — will skip merge.")
        return pd.DataFrame()

    traits = pd.read_csv(tsv_path, sep="\t", encoding="utf-8", low_memory=False)

    # identify genus column
    genus_col = None
    for c in traits.columns:
        if c.strip().lower() in {"genus", "género", "genero", "genus_clean", "genusname"}:
            genus_col = c
            break

    if genus_col is None and "GENUS" in traits.columns:
        genus_col = "GENUS"

    if genus_col is None:
        raise ValueError("Traits DB does not contain a recognizable GENUS column.")

    # standardize join key
    traits = traits.copy()
    traits["GENUS"] = traits[genus_col].astype(str)
    traits["GENUS_clean"] = traits["GENUS"].str.strip()

    # drop empty keys and deduplicate key
    traits = traits[traits["GENUS_clean"] != ""].copy()
    traits = traits.drop_duplicates(subset=["GENUS_clean"], keep="first").copy()

    print(f"🍄 Traits loaded: {len(traits)} unique GENUS_clean")
    return traits

def drop_mangled_duplicates_by_base(df: pd.DataFrame, keep_last_for: set[str]) -> pd.DataFrame:
    """
    If pandas created duplicated columns like primary_lifestyle and primary_lifestyle.1,
    drop earlier occurrences and keep the last for those in keep_last_for.
    For all other duplicated bases, keep the first.
    Then rename kept columns back to base name (remove .1).
    """
    cols = list(df.columns)

    def base(c):
        c = str(c).strip()
        return re.sub(r"\.\d+$", "", c)

    bases = [base(c) for c in cols]

    # mapping base -> positions
    pos = {}
    for i, b in enumerate(bases):
        pos.setdefault(b, []).append(i)

    keep_idx = set()
    removed_cols = []

    for b, idxs in pos.items():
        if len(idxs) == 1:
            keep_idx.add(idxs[0])
        else:
            if b in keep_last_for:
                keep_idx.add(idxs[-1])  # keep last
                removed_cols.extend([cols[i] for i in idxs[:-1]])
            else:
                keep_idx.add(idxs[0])   # keep first
                removed_cols.extend([cols[i] for i in idxs[1:]])

    if removed_cols:
        print(f"🧼 Removed duplicated/mangled columns: {removed_cols}")

    df = df.iloc[:, [i for i in range(len(cols)) if i in keep_idx]].copy()

    # rename mangled names to base
    rename_map = {c: base(c) for c in df.columns}
    df = df.rename(columns=rename_map)

    return df

def order_final_columns(df: pd.DataFrame, site_cols_final: list, traits_df: pd.DataFrame) -> pd.DataFrame:
    """
    Final order:
      fixed cols (BASE_FIXED_COLS)
      + site columns
      + traits columns (ALL columns from traits, except join key GENUS_clean)
    We also ensure lifestyles appear at the end (inside traits block), not in fixed.
    """
    # fixed (keep only those that exist)
    fixed = [c for c in BASE_FIXED_COLS if c in df.columns]

    # traits columns to append (as they appear in traits file)
    traits_cols = []
    if traits_df is not None and not traits_df.empty:
        traits_cols = [c for c in traits_df.columns if c != "GENUS_clean"]

    # remove any duplicates already present in fixed/sites
    traits_cols = [c for c in traits_cols if c not in fixed and c not in site_cols_final]

    final_cols = fixed + site_cols_final + traits_cols

    # keep any unexpected columns at the very end (just in case)
    extras = [c for c in df.columns if c not in final_cols]
    if extras:
        print(f"ℹ️ Extra columns appended at end (not in template): {extras}")
        final_cols += extras

    return df[final_cols].copy()

# ===================== MAIN =====================

def main():
    print("=" * 70)
    print("🚀 ATLASMX ITS — DB BUILDER V9.1 FULL PIPELINE")
    print("=" * 70)
    print('SUBSET BY TAXONOMY + PRUNE (OTUs & SITES) + DROP INPUT LIFESTYLES + MERGE FULL FUNGAL TRAITS')
    print("=" * 70)
    print()

    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input not found: {INPUT_FILE}")

    df = pd.read_csv(INPUT_FILE)
    print(f"📦 Loaded input: {len(df)} rows, {len(df.columns)} cols")

    # 1) Subset by taxonomy
    subset = apply_taxonomic_filter(df, TAXONOMIC_RANK, TAXONOMIC_FILTER, TAX_FILTER_MODE)

    # 2) Detect sites (from ORIGINAL subset structure)
    site_cols = detect_site_columns(subset)
    print(f"🔎 Site columns detected: {len(site_cols)}")

    # 3) Ensure numeric sites + prune zero-OTUs
    if site_cols:
        subset[site_cols] = subset[site_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

        if DROP_ZERO_OTUS:
            before = len(subset)
            subset = subset[subset[site_cols].sum(axis=1) > 0].copy()
            after = len(subset)
            print(f"🗜️ Pruned zero-OTUs by all sites: {before} → {after} (removed {before-after})")
    else:
        print("⚠️ No site columns detected. OTU pruning skipped.")

    # 4) Drop input lifestyles BEFORE merge (avoid duplicates)
    present_drop = [c for c in subset.columns if c in DROP_INPUT_LIFESTYLES]
    if present_drop:
        subset = subset.drop(columns=present_drop, errors="ignore")
        print(f"🧹 Dropped INPUT lifestyle cols (to avoid duplicates): {present_drop}")
    else:
        print("ℹ️ No input lifestyle columns found to drop.")

    # 5) Rename site columns (Sitio_AtlasEXP01 -> AtlasEXP01)
    site_cols_final = site_cols.copy()
    if REMOVE_SITE_PREFIX and site_cols_final:
        renamed = rename_sites(site_cols_final, REMOVE_SITE_PREFIX)
        rename_map = {old: new for old, new in zip(site_cols_final, renamed) if old != new}
        if rename_map:
            subset = subset.rename(columns=rename_map)
            print(f"🔁 Renamed site columns (removed prefix '{REMOVE_SITE_PREFIX}'): {len(rename_map)}")
        site_cols_final = renamed

    # 5.5) NEW STEP: Prune zero-SITES (columns) AFTER OTU prune and AFTER rename
    if DROP_ZERO_SITES and site_cols_final:
        # Make sure they are numeric (defensive: rename doesn't change dtype but keep safe)
        subset[site_cols_final] = subset[site_cols_final].apply(pd.to_numeric, errors="coerce").fillna(0)

        site_sums = subset[site_cols_final].sum(axis=0)
        zero_sites = site_sums[site_sums <= 0].index.tolist()

        if zero_sites:
            subset = subset.drop(columns=zero_sites, errors="ignore")
            # update list of site columns
            site_cols_final = [c for c in site_cols_final if c not in set(zero_sites)]
            print(f"🧯 Pruned zero-SITE columns: removed {len(zero_sites)} sites with sum=0")
        else:
            print("✅ No zero-SITE columns to prune (all sites have >0 total reads).")

        print(f"📍 Remaining site columns: {len(site_cols_final)}")
    else:
        print("ℹ️ Zero-SITE pruning skipped (disabled or no site columns).")

    # 6) Rename Total_Abundance -> Total_Abundance_<suffix>
    if "Total_Abundance" in subset.columns:
        new_total = f"Total_Abundance_{SUFFIX_TOTAL_ABUNDANCE}"
        subset = subset.rename(columns={"Total_Abundance": new_total})
        print(f"🔁 Renamed Total_Abundance → {new_total}")

    # 7) Load traits and merge (FULL)
    traits = load_traits(FUNGAL_TRAITS_DB)
    if not traits.empty and "genus" in subset.columns:
        subset["GENUS_clean"] = subset["genus"].astype(str).str.strip()
        merged = subset.merge(traits, on="GENUS_clean", how="left", suffixes=("", "_traits"))
        merged = merged.drop(columns=["GENUS_clean"], errors="ignore")

        # 8) If any mangled duplicates exist (primary_lifestyle / Secondary_lifestyle), keep LAST
        merged = drop_mangled_duplicates_by_base(
            merged,
            keep_last_for={"primary_lifestyle", "Secondary_lifestyle"}
        )

        subset = merged
        print(f"🧩 Traits merged. Final cols now: {len(subset.columns)}")
    else:
        print("⚠️ Traits merge skipped (missing traits or missing 'genus' column).")

    # 9) Final ordering: fixed + sites + traits
    subset = order_final_columns(subset, site_cols_final, traits)

    # 10) Save
    out_csv = OUTPUT_DIR / f"{OUTPUT_NAME}.csv"
    subset.to_csv(out_csv, index=False, encoding="utf-8")
    print(f"💾 Saved subset CSV: {out_csv.resolve()}")

    # Optional FASTA
    if "SPPN" in subset.columns and "sequence" in subset.columns:
        out_fasta = OUTPUT_DIR / f"{OUTPUT_NAME}.fasta"
        with open(out_fasta, "w", encoding="utf-8") as f:
            n = 0
            for _, r in subset.iterrows():
                sppn = str(r["SPPN"]).strip()
                seq = str(r["sequence"]).strip()
                if sppn and seq and sppn.lower() != "nan" and seq.lower() != "nan":
                    f.write(f">{sppn}\n{seq}\n")
                    n += 1
        print(f"🧬 Saved FASTA: {out_fasta.resolve()} ({n} seqs)")

    # Sanity checks
    for col in ["primary_lifestyle", "Secondary_lifestyle"]:
        if col in subset.columns:
            print(f"✅ Column present once: {col}")
        else:
            print(f"ℹ️ Column not present: {col} (maybe traits missing)")

    print("\n📊 DONE")
    print(f"   • Rows: {len(subset)}")
    print(f"   • Cols: {len(subset.columns)}")

if __name__ == "__main__":
    main()