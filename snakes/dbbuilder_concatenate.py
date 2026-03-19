# ==================== ATLASMX ITS — DB BUILDER V9.1 FULL PIPELINE (PART 1 + PART 2 + FOR_R EXPORT) ====================
#
# PART 1 (Concatenation):
# 1) Abundance table is MASTER (only OTUs present there are kept).
# 2) Join taxonomy (SINTAX) + sequences (FASTA).
# 3) Save:
#    - concatenated_tables/<PREFIX>_concatenated.csv
#    - concatenated_tables/Fungi_concatenated/<PREFIX>_fungi.csv
#    - concatenated_tables/Non_annotated/<PREFIX>_non_annotated.csv
#    - concatenated_tables/AllEuk_concatenated/<PREFIX>_all_eukaryotes.csv   (UNFILTERED full table)

import pandas as pd
from pathlib import Path
from collections import defaultdict
import time
import re
import sys

# ==================== Check arguments ====================

try:
    INPUT_DIR_name = sys.argv[1]
    OUT_PREFIX_name = sys.argv[2]
    DATASET_name = sys.argv[3]
except IndexError:
    print("""
ATLASMX ITS — DB BUILDER V9.1

Usage:
dbbuilder_concatenation.py <input_directory> <output_prefix> <dataset>

Input directory:
    string indicating the path of the input data, for example "~/atlas/dbbuilder/to_build/input"
Output prefix:
    string indicating the path for the DBBuilder output, for example "~/atlas/dbbuilder/db_test"
Dataset:
    string indicating the prefix for the dataset to analyze, for example "ATLASMXC"

Example:
python snakes/dbbuilder_concatenate.py ~/atlas/dbbuilder/to_build/input/ ~/atlas/dbbuilder/db_test ATLASMXC
    """)
    sys.exit()

# ====================
# CONFIGURATION
# ====================         

# Part 1 I/O
INPUT_DIR = Path(INPUT_DIR_name)
OUT_PREFIX = Path(OUT_PREFIX_name) 
OUTPUT_DIR = OUT_PREFIX / "concatenated_tables"
FUNGI_DIR = OUTPUT_DIR / "Fungi_concatenated"
NON_ANNOTATED_DIR = OUTPUT_DIR / "Non_annotated"
ALL_EUK_DIR = OUTPUT_DIR / "AllEuk_concatenated"
DATASET = DATASET_name

# Create required directories
for directory in [OUTPUT_DIR, FUNGI_DIR, NON_ANNOTATED_DIR, ALL_EUK_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

# =========================================================
# ==================== PART 1 FUNCTIONS ====================
# =========================================================

def read_fasta_sequences(fasta_path: Path) -> dict:
    """Read FASTA file and return a dict {OTU: sequence}."""
    print(f"  Reading FASTA file: {fasta_path.name}")
    sequences = {}
    current_otu = None
    current_sequence = []

    with open(fasta_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_otu is not None:
                    sequences[current_otu] = "".join(current_sequence)
                current_otu = line[1:]
                current_sequence = []
            else:
                current_sequence.append(line)

        if current_otu is not None:
            sequences[current_otu] = "".join(current_sequence)

    print(f"    Sequences loaded: {len(sequences)}")
    return sequences


def read_taxonomy_sintax(taxonomy_path: Path) -> pd.DataFrame:
    """Read SINTAX taxonomy file into a DataFrame indexed by OTU with 'sintax_taxonomy' column."""
    print(f"  Reading SINTAX taxonomy: {taxonomy_path.name}")
    recs = []
    with open(taxonomy_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                recs.append({"OTU": parts[0], "sintax_taxonomy": parts[1]})
    df = pd.DataFrame(recs)
    if not df.empty:
        df = df.set_index("OTU")
    print(f"    Taxonomy records loaded: {len(df)}")
    return df


def process_single_dataset(prefix: str, files: dict) -> pd.DataFrame:
    """Build concatenated table (master=abundance) and sort by Total_Abundance."""
    print(f"\n=== PROCESSING DATASET: {prefix} ===")

    print("  Reading abundance table (master table)...")
    abundance_df = pd.read_csv(files["abundance"], dtype={"OTU": str}).set_index("OTU")
    print(f"    OTUs in abundance table: {len(abundance_df)}")

    taxonomy_df = read_taxonomy_sintax(files["taxonomy"])
    if not taxonomy_df.empty:
        merged_df = abundance_df.join(taxonomy_df, how="left")
    else:
        merged_df = abundance_df.copy()
        merged_df["sintax_taxonomy"] = ""

    sequences_dict = read_fasta_sequences(files["sequences"])
    merged_df["sequence"] = merged_df.index.map(sequences_dict)

    merged_df = merged_df.reset_index()
    merged_df["OTU_XX"] = [f"{prefix}_OTU{str(i+1).zfill(5)}" for i in range(len(merged_df))]
    merged_df["Original_OTUID"] = merged_df["OTU"]

    print("  Calculating and sorting by total abundance...")
    special = ["OTU", "sintax_taxonomy", "sequence", "OTU_XX", "Original_OTUID"]
    site_columns = [c for c in merged_df.columns if c not in special]

    merged_df["Total_Abundance"] = merged_df[site_columns].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1)
    merged_df = merged_df.sort_values("Total_Abundance", ascending=False)

    final_cols = ["OTU_XX", "Original_OTUID", "sintax_taxonomy", "sequence"] + site_columns
    final_df = merged_df[final_cols].copy()

    print(f"✅ Dataset '{prefix}' processed. Rows: {len(final_df)}")
    return final_df


def filter_and_save_datasets(concatenated_df: pd.DataFrame, prefix: str):
    """Save fungi/non-fungi subsets and the unfiltered all-eukaryote table."""
    print(f"=== FILTERING & SAVING DATASET: {prefix} ===")

    fungi_mask = concatenated_df["sintax_taxonomy"].astype(str).str.contains("d:Fungi", na=False)
    fungi_df = concatenated_df[fungi_mask].copy()
    non_df = concatenated_df[~fungi_mask].copy()

    fungi_path = FUNGI_DIR / f"{prefix}_fungi.csv"
    non_path = NON_ANNOTATED_DIR / f"{prefix}_non_annotated.csv"
    all_euk_path = ALL_EUK_DIR / f"{prefix}_all_eukaryotes.csv"

    fungi_df.to_csv(fungi_path, index=False)
    non_df.to_csv(non_path, index=False)
    concatenated_df.to_csv(all_euk_path, index=False)

    stats = {
        "prefix": prefix,
        "total_otus": len(concatenated_df),
        "fungi_otus": len(fungi_df),
        "non_annotated_otus": len(non_df),
        "fungi_percentage": (len(fungi_df) / len(concatenated_df) * 100) if len(concatenated_df) else 0.0,
    }

    print(f"✅ {prefix}: {stats['fungi_otus']} fungi, {stats['non_annotated_otus']} non-fungi ({stats['fungi_percentage']:.1f}%)")
    print(f"   💾 Saved ALL EUK file (unfiltered): {all_euk_path}")
    return stats


def generate_final_report(all_stats: list, processing_time: float, dataset: str):
    """Write final report log for Part 1."""
    print("\n" + "=" * 60)
    print("📊 FINAL PROCESSING REPORT (PART 1)")
    print("=" * 60)

    log_path = OUTPUT_DIR / f"{dataset}_processing_log.txt"
    with open(log_path, "w", encoding="utf-8") as log_file:
        log_file.write("PROCESSING REPORT - METABARCODING (PART 1)\n")
        log_file.write("Logic: Abundance master + sort by total abundance\n")
        log_file.write("=" * 60 + "\n\n")
        for s in all_stats:
            log_file.write(f"Dataset: {s['prefix']}\n")
            log_file.write(f"  Total OTUs: {s['total_otus']}\n")
            log_file.write(f"  Fungi OTUs: {s['fungi_otus']} ({s['fungi_percentage']:.1f}%)\n")
            log_file.write(f"  Non-fungi OTUs: {s['non_annotated_otus']}\n\n")
        log_file.write(f"Total processing time: {processing_time:.2f} seconds\n")

    total_otus = sum(s["total_otus"] for s in all_stats)
    total_fungi = sum(s["fungi_otus"] for s in all_stats)
    print(f"📁 Processed datasets: {len(all_stats)}")
    print(f"🧬 Total OTUs processed: {total_otus}")
    print(f"🍄 Total Fungi OTUs: {total_fungi} ({(total_fungi/total_otus*100 if total_otus else 0):.1f}%)")
    print(f"⏱️ Total time: {processing_time:.2f} seconds")
    print(f"📋 Detailed log saved at: {log_path}")

# =========================================================
# ==================== MAIN (PART 1) ======================
# =========================================================

def main():
    print("=" * 70)
    print("🚀 ATLASMX ITS — DB BUILDER V9.1")
    print("=" * 70)
    print("PART 1: Concatenate + save Fungi / Non-fungi / AllEuk")
    print("=" * 70)
    print()

    # ---------------- PART 1 ----------------
    print("\n🧩 PART 1 — STARTING METABARCODING CONCATENATION")
    start_time = time.time()

    datasets =  {DATASET:{'abundance': INPUT_DIR / f'{DATASET}_otu_abundance_table.csv',
                         'sequences': INPUT_DIR / f'{DATASET}_sequences.fasta',
                         'taxonomy': INPUT_DIR / f'{DATASET}_taxonomy.sintax.txt'}}
    if not datasets:
        print("❌ No complete datasets found in 'input' folder")
        return

    all_stats = []
    for prefix, files in datasets.items():
        try:
            concatenated_df = process_single_dataset(prefix, files)

            concat_output_path = OUTPUT_DIR / f"{prefix}_concatenated.csv"
            concatenated_df.to_csv(concat_output_path, index=False)
            print(f"✅ Concatenated saved: {concat_output_path}")

            stats = filter_and_save_datasets(concatenated_df, prefix)
            all_stats.append(stats)

        except Exception as e:
            print(f"❌ Error processing dataset {prefix}: {e}")
            continue

    processing_time = time.time() - start_time
    generate_final_report(all_stats, processing_time, DATASET)
    

if __name__ == "__main__":
    main()