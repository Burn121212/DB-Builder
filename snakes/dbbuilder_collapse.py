# ==================== ATLASMX ITS — DB BUILDER V9.1 FULL PIPELINE (PART 2) ====================
# 
# PART 2 (Taxonomic collapse + final DB):
# 4) Collapse according to COLLAPSE_STRATEGY:
#    - "species_only" -> collapse only when species p >= P_VALUE_THRESHOLD
#    - "genus"        -> collapse only when genus p >= P_VALUE_THRESHOLD
#    - "all"          -> recursive lowest-rank collapse using the deepest available rank
#                        with p >= P_VALUE_THRESHOLD among:
#                        species -> genus -> family -> order -> class -> phylum -> domain
# 5) Representative is chosen ONLY among rows with COMPLETE taxonomy (d,p,c,o,f,g,s);
#    among those, the longest sequence is kept.
# 6) Sum abundances across sites.
# 7) Compute Total_Abundance, sort, and build binomial species names.
# 8) SANITIZE taxonomy columns by p >= SPPN_P_THRESHOLD
#    (fix SINTAX cases where names remain despite low confidence).
# 9) Add fungal lifestyles (Fungal_Traits_DB) after sanitization.
# 10) Assign SPPN (deepest rank with p >= SPPN_P_THRESHOLD), prune zeros.
# 11) Save final DB CSVs with strategy and p-value in the filename:
#     - FINAL_DB/Final_Database_<strategy>_all_eukaryotes_p{P_VALUE_THRESHOLD}.csv   (when MODE=all_eukaryotes)
#     - FINAL_DB/Final_Database_<strategy>_fungi_p{P_VALUE_THRESHOLD}.csv             (subset from all-euk OR main when MODE=fungi)
#
#
# NOTE:
# - MODE controls which tables are used for collapse in Part 2.
# - Regardless of MODE, if MODE="all_eukaryotes" we also export a fungi-only final DB + its For_R package.
# - This version preserves the original conservative logic for representative selection:
#   a group is collapsed only if at least one row in that group has COMPLETE taxonomy.
#   Otherwise, the group is left uncollapsed.

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
ATLASMX ITS — DB BUILDER V9.1 COLLAPSE

Usage:
dbbuilder_concatenation.py <output_prefix> <collapse_mode> <collapse_strategy> <p-value_threshold>

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
python snakes/dbbuilder_concatenate.py ~/atlas/dbbuilder/db_test all_eukaryotes species_only 1.0
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
for directory in [OUTPUT_DIR, FUNGI_DIR, NON_ANNOTATED_DIR, ALL_EUK_DIR, FINAL_DB_DIR]:
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


def build_collapse_log_filename(suffix: str) -> str:
    """
    Build collapse log filename using:
    strategy + suffix + collapse threshold
    """
    return f"Collapse_Log_{COLLAPSE_STRATEGY}_{suffix}_p{P_VALUE_THRESHOLD}.txt"

# =========================================================
# ==================== PART 2 FUNCTIONS ====================
# =========================================================

def clean_taxonomy_string(tax_str):
    """Remove artifacts such as '(Fungi)' inserted inside taxon names."""
    if not isinstance(tax_str, str):
        return tax_str
    return re.sub(r"([kpcofgs]):([^(),]+)\([^)]+\)(?=\(\d+\.\d+\))", r"\1:\2", tax_str)


def extract_taxonomy_data(tax_str):
    """
    Parse a SINTAX string into a dict:
    {rank_prefix: {'taxon': ..., 'p_value': ...}}
    """
    if not isinstance(tax_str, str):
        return {}
    clean_str = clean_taxonomy_string(tax_str)
    tax_data = {}
    for prefix, taxon, p_value in re.findall(r"([dkpcofgs]):([^(]+)\(([\d.]+)\)", clean_str):
        tax_data[prefix] = {"taxon": taxon.strip(), "p_value": float(p_value)}
    return tax_data


def build_collapse_key_species_only(tax_data, p_thresh):
    """Build collapse key only if species exists and species p-value >= threshold."""
    if not tax_data or "s" not in tax_data:
        return None
    if tax_data["s"]["p_value"] < p_thresh:
        return None

    hierarchy = ["d", "p", "c", "o", "f", "g", "s"]
    parts = []
    for level in hierarchy:
        if level in tax_data:
            parts.append(f"{level}:{tax_data[level]['taxon']}")
        if level == "s":
            break
    return "_".join(parts) if parts else None


def build_collapse_key_genus(tax_data, p_thresh):
    """Build collapse key only if genus exists and genus p-value >= threshold."""
    if not tax_data or "g" not in tax_data:
        return None
    if tax_data["g"]["p_value"] < p_thresh:
        return None

    hierarchy = ["d", "p", "c", "o", "f", "g"]
    parts = []
    for level in hierarchy:
        if level in tax_data:
            parts.append(f"{level}:{tax_data[level]['taxon']}")
        if level == "g":
            break
    return "_".join(parts) if parts else None


def build_collapse_key_lowest_rank(tax_data, p_thresh):
    """
    Build collapse key using the deepest available confident rank among:
    species -> genus -> family -> order -> class -> phylum -> domain

    The resulting key includes all ranks available up to the chosen level.
    """
    if not tax_data:
        return None

    ranked_levels = ["s", "g", "f", "o", "c", "p", "d"]
    chosen = None

    for level in ranked_levels:
        if level in tax_data and tax_data[level]["p_value"] >= p_thresh:
            chosen = level
            break

    if chosen is None:
        return None

    hierarchy = ["d", "p", "c", "o", "f", "g", "s"]
    parts = []
    for level in hierarchy:
        if level in tax_data:
            parts.append(f"{level}:{tax_data[level]['taxon']}")
        if level == chosen:
            break

    return "_".join(parts) if parts else None


def build_collapse_key(tax_data, strategy, p_thresh):
    """Dispatch collapse key generation according to the selected strategy."""
    if strategy == "species_only":
        return build_collapse_key_species_only(tax_data, p_thresh)
    elif strategy == "genus":
        return build_collapse_key_genus(tax_data, p_thresh)
    elif strategy == "all":
        return build_collapse_key_lowest_rank(tax_data, p_thresh)
    else:
        raise ValueError(f"Unsupported COLLAPSE_STRATEGY: {strategy}")


def parse_tax_columns_from_sintax(sintax_str):
    """Build explicit taxonomy columns and p-value columns from a SINTAX string."""
    tax = extract_taxonomy_data(sintax_str)
    ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    prefs = ["d", "p", "c", "o", "f", "g", "s"]
    out = {}
    for r, p in zip(ranks, prefs):
        out[r] = tax.get(p, {}).get("taxon", "")
        out[f"{r}_pvalue"] = tax.get(p, {}).get("p_value", 0.0)
    return out


def has_complete_taxonomy(sintax_str):
    """Return True only if d,p,c,o,f,g,s all exist and are non-empty."""
    tax = extract_taxonomy_data(sintax_str)
    required = ["d", "p", "c", "o", "f", "g", "s"]
    for k in required:
        if k not in tax or not str(tax[k].get("taxon", "")).strip():
            return False
    return True


def get_confident_rank_counts(tx_series, p_thresh):
    """
    Return counts of rows having a confident assignment at each rank.
    Used for strategy-aware reporting in Part 2.
    """
    return {
        "species": int(tx_series.apply(lambda x: ("s" in x and x["s"]["p_value"] >= p_thresh) if x else False).sum()),
        "genus": int(tx_series.apply(lambda x: ("g" in x and x["g"]["p_value"] >= p_thresh) if x else False).sum()),
        "family": int(tx_series.apply(lambda x: ("f" in x and x["f"]["p_value"] >= p_thresh) if x else False).sum()),
        "order": int(tx_series.apply(lambda x: ("o" in x and x["o"]["p_value"] >= p_thresh) if x else False).sum()),
        "class": int(tx_series.apply(lambda x: ("c" in x and x["c"]["p_value"] >= p_thresh) if x else False).sum()),
        "phylum": int(tx_series.apply(lambda x: ("p" in x and x["p"]["p_value"] >= p_thresh) if x else False).sum()) if False else 0
    }


def get_deepest_confident_rank_label(tax_data, p_thresh):
    """
    Return the deepest confident rank label among:
    species > genus > family > order > class > phylum > domain
    """
    if not tax_data:
        return None

    if "s" in tax_data and tax_data["s"]["p_value"] >= p_thresh:
        return "species"
    if "g" in tax_data and tax_data["g"]["p_value"] >= p_thresh:
        return "genus"
    if "f" in tax_data and tax_data["f"]["p_value"] >= p_thresh:
        return "family"
    if "o" in tax_data and tax_data["o"]["p_value"] >= p_thresh:
        return "order"
    if "c" in tax_data and tax_data["c"]["p_value"] >= p_thresh:
        return "class"
    if "p" in tax_data and tax_data["p"]["p_value"] >= p_thresh:
        return "phylum"
    if "d" in tax_data and tax_data["d"]["p_value"] >= p_thresh:
        return "domain"
    return None


def load_datasets_with_stats_part2():
    """Load Part 2 input tables according to MODE and print strategy-aware summary statistics."""
    if MODE == "fungi":
        input_dir = FUNGI_DIR
        pattern = "*_fungi.csv"
    else:
        input_dir = ALL_EUK_DIR
        pattern = "*_all_eukaryotes.csv"

    print("\n📁 PART 2 — LOADING FILES...")
    files = list(input_dir.glob(pattern))
    if not files:
        print(f"❌ No files found in: {input_dir} with pattern: {pattern}")
        return pd.DataFrame(), [], input_dir, pattern

    all_dfs = []
    site_union = set()
    subset_stats = []

    fixed = {"OTU_XX", "Original_OTUID", "sintax_taxonomy", "sequence"}

    for csv_file in files:
        df = pd.read_csv(csv_file)
        subset = csv_file.stem.replace("_fungi", "").replace("_all_eukaryotes", "")
        site_cols = [c for c in df.columns if c not in fixed]

        print(f"   🧹 Cleaning taxonomy strings for {subset}...")
        df["sintax_taxonomy"] = df["sintax_taxonomy"].apply(clean_taxonomy_string)

        tx = df["sintax_taxonomy"].apply(extract_taxonomy_data)
        complete_tax = int(df["sintax_taxonomy"].apply(has_complete_taxonomy).sum())

        with_species = int(tx.apply(lambda x: ("s" in x) if x else False).sum())
        with_genus = int(tx.apply(lambda x: ("g" in x) if x else False).sum())

        confident_species = int(tx.apply(lambda x: ("s" in x and x["s"]["p_value"] >= P_VALUE_THRESHOLD) if x else False).sum())
        confident_genus = int(tx.apply(lambda x: ("g" in x and x["g"]["p_value"] >= P_VALUE_THRESHOLD) if x else False).sum())

        deepest_counts = {
            "species": 0,
            "genus": 0,
            "family": 0,
            "order": 0,
            "class": 0,
            "phylum": 0,
            "domain": 0,
            "none": 0
        }

        for tax_data in tx:
            label = get_deepest_confident_rank_label(tax_data, P_VALUE_THRESHOLD)
            if label is None:
                deepest_counts["none"] += 1
            else:
                deepest_counts[label] += 1

        subset_stats.append({
            "subset": subset,
            "otus": len(df),
            "sites": len(site_cols),
            "with_species": with_species,
            "with_genus": with_genus,
            "confident_species": confident_species,
            "confident_genus": confident_genus,
            "complete_tax": complete_tax,
            "species_%": (with_species / len(df) * 100) if len(df) else 0.0,
            "genus_%": (with_genus / len(df) * 100) if len(df) else 0.0,
            "complete_tax_%": (complete_tax / len(df) * 100) if len(df) else 0.0,
            "deepest_counts": deepest_counts
        })

        site_union.update(site_cols)
        all_dfs.append((subset, df))

    print("\n📊 PART 2 — SUBSET STATISTICS")
    print("-" * 120)

    if COLLAPSE_STRATEGY == "species_only":
        print(f"{'Subset':<24}{'OTUs':>8}{'Sites':>8}{'with s':>10}{'s p≥p':>10}{'complete':>10}{'% s':>8}{'% complete':>12}")
        print("-" * 120)
        for s in subset_stats:
            print(
                f"{s['subset']:<24}{s['otus']:>8}{s['sites']:>8}{s['with_species']:>10}"
                f"{s['confident_species']:>10}{s['complete_tax']:>10}{s['species_%']:>8.1f}{s['complete_tax_%']:>12.1f}"
            )

    elif COLLAPSE_STRATEGY == "genus":
        print(f"{'Subset':<24}{'OTUs':>8}{'Sites':>8}{'with g':>10}{'g p≥p':>10}{'complete':>10}{'% g':>8}{'% complete':>12}")
        print("-" * 120)
        for s in subset_stats:
            print(
                f"{s['subset']:<24}{s['otus']:>8}{s['sites']:>8}{s['with_genus']:>10}"
                f"{s['confident_genus']:>10}{s['complete_tax']:>10}{s['genus_%']:>8.1f}{s['complete_tax_%']:>12.1f}"
            )

    else:
        print(
            f"{'Subset':<24}{'OTUs':>8}{'Sites':>8}{'species':>10}{'genus':>10}{'family':>10}"
            f"{'order':>10}{'class':>10}{'phylum':>10}{'domain':>10}{'none':>10}"
        )
        print("-" * 120)
        for s in subset_stats:
            dc = s["deepest_counts"]
            print(
                f"{s['subset']:<24}{s['otus']:>8}{s['sites']:>8}"
                f"{dc['species']:>10}{dc['genus']:>10}{dc['family']:>10}"
                f"{dc['order']:>10}{dc['class']:>10}{dc['phylum']:>10}{dc['domain']:>10}{dc['none']:>10}"
            )

    site_cols_all = sorted(site_union)

    unified = []
    for subset, df in all_dfs:
        missing = set(site_cols_all) - set([c for c in df.columns if c not in fixed])
        if missing:
            df_missing = pd.DataFrame(0, index=df.index, columns=sorted(missing))
            df = pd.concat([df, df_missing], axis=1)

        df = df[["OTU_XX", "Original_OTUID", "sintax_taxonomy", "sequence"] + site_cols_all]
        unified.append(df)
        print(f"   ✅ {subset}: {len(df)} OTUs")

    combined = pd.concat(unified, ignore_index=True, copy=False)
    print("\n📈 PART 2 — COMBINED SUMMARY")
    print(f"   • Total OTUs: {len(combined)}")
    print(f"   • Total sites: {len(site_cols_all)}")
    print(f"   • MODE: {MODE} | Input: {input_dir} | Pattern: {pattern}")
    print(f"   • COLLAPSE_STRATEGY: {COLLAPSE_STRATEGY}")
    print(f"   • P_VALUE_THRESHOLD: {P_VALUE_THRESHOLD}")
    print(f"   • SPPN_P_THRESHOLD: {SPPN_P_THRESHOLD}")

    return combined, site_cols_all, input_dir, pattern


def create_binomial_species_inplace(df):
    """
    Convert species labels to Genus_species when both genus and species are available.
    This is done in-place after collapse and before SPPN assignment.
    """
    created = 0
    no_genus = 0
    no_species = 0

    for idx, row in df.iterrows():
        genus = row.get("genus", "")
        species = row.get("species", "")

        if genus and species:
            df.at[idx, "species"] = f"{genus}_{species}"
            created += 1
        elif species and not genus:
            no_genus += 1
        elif genus and not species:
            no_species += 1

    print(f"   • Binomials created: {created}")
    if no_genus or no_species:
        print(f"   • Warning: no genus for species: {no_genus} | no species for genus: {no_species}")


def sanitize_taxonomy_columns_by_threshold(df: pd.DataFrame, threshold: float, cascade_blank: bool = True) -> pd.DataFrame:
    """
    Blank taxon names when p-value < threshold.

    If cascade_blank=True, once a higher rank is invalid,
    all deeper ranks are also blanked.
    """
    print(f"\n🧽 SANITIZING TAXONOMY COLUMNS BY p≥{threshold:.2f} (cascade_blank={cascade_blank})...")
    ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]

    for r in ranks:
        pcol = f"{r}_pvalue"
        if pcol in df.columns:
            df[pcol] = pd.to_numeric(df[pcol], errors="coerce").fillna(0.0)
        else:
            df[pcol] = 0.0
        if r not in df.columns:
            df[r] = ""

    blank_counts = {r: 0 for r in ranks}

    if not cascade_blank:
        for r in ranks:
            mask = df[f"{r}_pvalue"] < threshold
            blank_counts[r] = int((mask & (df[r].astype(str).str.strip() != "")).sum())
            df.loc[mask, r] = ""
    else:
        invalid_up_to = pd.Series(False, index=df.index)
        for r in ranks:
            invalid_up_to = invalid_up_to | (df[f"{r}_pvalue"] < threshold)
            mask = invalid_up_to
            blank_counts[r] = int((mask & (df[r].astype(str).str.strip() != "")).sum())
            df.loc[mask, r] = ""

    total_blank = sum(blank_counts.values())
    if total_blank:
        print("   • Blanked names per rank:")
        for r in ranks:
            if blank_counts[r]:
                print(f"     - {r}: {blank_counts[r]}")
    else:
        print("   • Nothing blanked (no low-p names present).")

    return df


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


def add_lifestyle_info(final_df, fungal_traits):
    """Add primary and secondary lifestyles by exact genus match (case-insensitive)."""
    print("\n🔄 ADDING LIFESTYLE INFORMATION...")
    final_df["primary_lifestyle"] = ""
    final_df["Secondary_lifestyle"] = ""

    if fungal_traits is None or fungal_traits.empty:
        print("   • No fungal traits available — leaving empty lifestyle columns.")
        return final_df

    lookup = {}
    for _, r in fungal_traits.iterrows():
        g = str(r["GENUS_clean"]).strip()
        if g:
            lookup[g.lower()] = (r["primary_lifestyle"], r["Secondary_lifestyle"])

    matches = 0
    for i, row in final_df.iterrows():
        genus = str(row.get("genus", "")).strip()
        if genus:
            key = genus.lower()
            if key in lookup:
                p, s = lookup[key]
                final_df.at[i, "primary_lifestyle"] = p
                final_df.at[i, "Secondary_lifestyle"] = s
                matches += 1

    print(f"   • Lifestyle filled for {matches} rows ({(matches/len(final_df)*100 if len(final_df) else 0):.1f}%)")
    return final_df


def collapse_otus_simple(combined_df, site_cols):
    """
    Collapse OTUs using COLLAPSE_STRATEGY:
      - species_only
      - genus
      - all (deepest available confident rank)

    The original conservative representative selection rule is preserved:
    representative is chosen ONLY among rows with COMPLETE taxonomy (d,p,c,o,f,g,s),
    and the longest sequence is selected.
    """
    print(f"\n🔬 APPLYING COLLAPSE STRATEGY: {COLLAPSE_STRATEGY} (p≥{P_VALUE_THRESHOLD})")

    tax_list = combined_df["sintax_taxonomy"].apply(extract_taxonomy_data)
    collapse_keys = tax_list.apply(lambda t: build_collapse_key(t, COLLAPSE_STRATEGY, P_VALUE_THRESHOLD))

    total = len(combined_df)
    collapsable = int(collapse_keys.notna().sum())
    non_collapsable = total - collapsable

    complete_in_collapsable = 0
    for idx in collapse_keys[collapse_keys.notna()].index:
        if has_complete_taxonomy(combined_df.loc[idx, "sintax_taxonomy"]):
            complete_in_collapsable += 1

    print("📊 OTU DISTRIBUTION:")
    print(f"   • Collapsable: {collapsable} ({(collapsable/total*100 if total else 0):.1f}%)")
    print(f"   • Complete among collapsable: {complete_in_collapsable} ({(complete_in_collapsable/collapsable*100 if collapsable else 0):.1f}%)")
    print(f"   • Individual: {non_collapsable} ({(non_collapsable/total*100 if total else 0):.1f}%)")

    result_rows = []
    logs = []

    group_index = {}
    for idx, key in collapse_keys.items():
        group_index.setdefault(key, []).append(idx)

    print(f"\n🔄 Processing {len([k for k in group_index if k is not None])} collapsable groups...")

    for key, idxs in group_index.items():
        grp = combined_df.loc[idxs]

        if key is None:
            for _, row in grp.iterrows():
                parsed = parse_tax_columns_from_sintax(row["sintax_taxonomy"])
                out = row.to_dict()
                out.update(parsed)
                result_rows.append(out)
            continue

        if len(grp) == 1:
            row = grp.iloc[0]
            parsed = parse_tax_columns_from_sintax(row["sintax_taxonomy"])
            out = row.to_dict()
            out.update(parsed)
            result_rows.append(out)
            continue

        grp = grp.copy()
        complete_mask = grp["sintax_taxonomy"].apply(has_complete_taxonomy)
        candidates = grp[complete_mask]

        if candidates.empty:
            for _, row in grp.iterrows():
                parsed = parse_tax_columns_from_sintax(row["sintax_taxonomy"])
                out = row.to_dict()
                out.update(parsed)
                result_rows.append(out)

            logs.append({
                "collapse_key": key,
                "kept_otu": None,
                "removed_count": 0,
                "removed_otus": [],
                "group_size": len(grp),
                "representative_seq_len": None,
                "representative_tax": None,
                "note": f"skip_collapse_incomplete_tax_{COLLAPSE_STRATEGY}",
            })
            continue

        seq_lengths = candidates["sequence"].astype(str).str.len().fillna(0)
        idx_keep = seq_lengths.idxmax()
        kept = candidates.loc[idx_keep].copy()

        summed = grp[site_cols].apply(pd.to_numeric, errors="coerce").fillna(0).sum()

        parsed_tax = parse_tax_columns_from_sintax(kept["sintax_taxonomy"])
        out = kept.to_dict()
        out.update(parsed_tax)
        for c, val in zip(site_cols, summed.values):
            out[c] = val
        result_rows.append(out)

        removed = [r for r in grp["OTU_XX"].tolist() if r != kept["OTU_XX"]]
        logs.append({
            "collapse_key": key,
            "kept_otu": kept["OTU_XX"],
            "removed_count": len(removed),
            "removed_otus": removed[:10],
            "group_size": len(grp),
            "representative_seq_len": int(seq_lengths.loc[idx_keep]),
            "representative_tax": kept["sintax_taxonomy"],
            "note": f"collapsed_with_complete_tax_{COLLAPSE_STRATEGY}",
        })

    final_df = pd.DataFrame(result_rows)
    print("\n✅ COLLAPSE COMPLETED")
    print(f"   • Initial OTUs: {len(combined_df)}")
    print(f"   • Final OTUs:   {len(final_df)}")
    print(f"   • Reduction: {len(combined_df) - len(final_df)} OTUs")

    final_df["sequence_lenght"] = final_df["sequence"].astype(str).str.len()
    return final_df, logs


def compute_total_abundance_and_sort(final_df, site_cols):
    """Compute Total_Abundance, sort descending, and build binomial species labels."""
    print("\n📊 COMPUTING TOTAL ABUNDANCE & SORTING...")
    site_data = final_df[site_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    final_df["Total_Abundance"] = site_data.sum(axis=1)
    final_df = final_df.sort_values("Total_Abundance", ascending=False).reset_index(drop=True)
    create_binomial_species_inplace(final_df)
    return final_df


def assign_sppn_unique(final_df):
    """
    Assign SPPN codes using the deepest rank whose p-value is >= SPPN_P_THRESHOLD.
    Priority:
    species > genus > family > order > class > phylum > domain
    """
    print("\n🏷️ ASSIGNING FINAL SPPN CODES (PER TAXON, p≥{:.2f})...".format(SPPN_P_THRESHOLD))
    if "SPPN" in final_df.columns:
        final_df = final_df.drop(columns=["SPPN"])

    sppn = []
    counter = {}
    rank_order = ["species", "genus", "family", "order", "class", "phylum", "domain"]

    for _, row in final_df.iterrows():
        base_taxon = None
        for rank in rank_order:
            name = str(row.get(rank, "") or "").strip()
            pval = row.get(f"{rank}_pvalue", 0.0)
            try:
                pval = float(pval)
            except Exception:
                pval = 0.0
            if name and pval >= SPPN_P_THRESHOLD:
                base_taxon = name
                break

        if not base_taxon:
            base_taxon = "Unknown"

        base = re.sub(r"[^\w]+", "_", base_taxon)
        counter[base] = counter.get(base, 0) + 1
        sppn.append(f"{base}_{counter[base]:04d}")

    final_df["SPPN"] = sppn
    return final_df


def prune_zero_abundance(final_df, site_cols):
    """Drop rows with Total_Abundance == 0."""
    print("\n🧹 PRUNING ZERO-ABUNDANCE OTUs...")
    site_data = final_df[site_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    final_df["Total_Abundance"] = site_data.sum(axis=1)
    zero_mask = final_df["Total_Abundance"] == 0
    n_zero = int(zero_mask.sum())

    if n_zero:
        print(f"   • Removing {n_zero} OTUs with Total_Abundance == 0")
        final_df = final_df.loc[~zero_mask].reset_index(drop=True)
    else:
        print("   • No zeros to prune.")

    print(f"   • Remaining OTUs: {len(final_df)}")
    return final_df


TAX_RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
PVALUE_COLS = [f"{r}_pvalue" for r in TAX_RANKS]

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
# ==================== MAIN (PART 2) ======================
# =========================================================

def main():
    print("=" * 70)
    print("🚀 ATLASMX ITS — DB BUILDER V9.1")
    print("=" * 70)
    print(f"PART 2: Collapse ({COLLAPSE_STRATEGY}) + final database(s) in FINAL_DB")
    print(f"MODE: {MODE} | strategy: {COLLAPSE_STRATEGY} | collapse p≥{P_VALUE_THRESHOLD} | sanitize p≥{SPPN_P_THRESHOLD}")
    print("=" * 70)
    print()    

    # ---------------- PART 2 ----------------
    print("\n🧬 PART 2 — STARTING TAXONOMIC COLLAPSE")
    fungal_traits = load_fungal_traits_db()

    combined, site_cols, part2_input_dir, part2_pattern = load_datasets_with_stats_part2()
    if combined.empty:
        return

    final_df, logs = collapse_otus_simple(combined, site_cols)
    final_df = compute_total_abundance_and_sort(final_df, site_cols)

    final_df = sanitize_taxonomy_columns_by_threshold(
        final_df,
        threshold=SPPN_P_THRESHOLD,
        cascade_blank=True
    )

    final_df = add_lifestyle_info(final_df, fungal_traits)
    final_df = assign_sppn_unique(final_df)
    final_df = prune_zero_abundance(final_df, site_cols)

    produced_csvs = save_outputs(final_df, logs, site_cols, part2_input_dir, part2_pattern)


    print("\n🎉 ALL DONE")
    print(f"   MODE: {MODE}")
    print(f"   COLLAPSE_STRATEGY: {COLLAPSE_STRATEGY}")
    print(f"   Final OTUs: {len(final_df)}")
    print(f"   FINAL_DB folder: {FINAL_DB_DIR.resolve()}")

if __name__ == "__main__":
    main()