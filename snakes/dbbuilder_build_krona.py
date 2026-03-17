# ============================================================
# BUILD KRONA
# ============================================================

from pathlib import Path
import pandas as pd
import sys

# ==================== Check arguments ====================

try:
    OUT_PREFIX_name = sys.argv[1]
    INPUT_File_name = sys.argv[2]
    OUTPUT_name = sys.argv[3]
    DATASET_name = sys.argv[4]
    INCLUDE_BLANK_COLUMN_in = sys.argv[5]
    
except IndexError:
    print("""
ATLASMX ITS — DB BUILDER V9.1
BUILD KRONA

Usage:
dbbuilder_build_krona.py <output_prefix> <input_file_name> <output_name> <dataset_name> <include_blank_column>

Output prefix:
    string indicating the path of the DBBuilder output, for example "~/atlas/dbbuilder/db_test"
Input file name:
    string for the input file, for example "FINAL_DB/Final_Database_species_only_fungi_p0.8.csv"
Output name:
    string for output  file name, for example "Krona_format_species_only_fungi_p0.8_zOTU.csv"
Dataset name:
    string to repeat in the first column, for example "fungi_zOTU_p0.8"
Include blank column (between amount and taxonomy):
    0    ->    no
    1    ->    yes

Example:
python snakes/dbbuilder_build_krona.py ~/atlas/dbbuilder/db_test FINAL_DB/Final_Database_species_only_fungi_p1.0.csv Krona_format_species_only_fungi_p1.0_zOTU.csv fungi_zOTU_p1.0 1
    """)
    sys.exit()


# =========================
# CONFIGURATION 
# =========================
OUT_PREFIX = Path(OUT_PREFIX_name) # latest version change
INPUT_FILE = OUT_PREFIX / INPUT_File_name # latest version change

# carpeta de salida
OUTPUT_DIR = OUT_PREFIX / "krona2" # latest version change

# archivo de salida
OUTPUT_FILE = OUTPUT_DIR / OUTPUT_name # latest version change

# nombre que quieres repetir en la primera columna
DATASET_NAME = DATASET_name

# si quieres incluir una columna vacía entre amount y taxonomía
INCLUDE_BLANK_COLUMN = bool(int(INCLUDE_BLANK_COLUMN_in))

def main():
    print("=" * 70)
    print("🚀 ATLASMX ITS — DB BUILDER V9.1 FULL PIPELINE")
    print("=" * 70)
    print('BUILD KRONA')
    print("=" * 70)
    print()

    # =========================
    # CREATE OUTPUT DIRECTORY
    # =========================
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Output directory ready: {OUTPUT_DIR.resolve()}")
    
    # =========================
    # LOAD INPUT
    # =========================
    df = pd.read_csv(INPUT_FILE, low_memory=False)
    
    print(f"Loaded input file: {INPUT_FILE}")
    print(f"Rows: {len(df):,}")
    print(f"Columns: {len(df.columns):,}")
    
    # =========================
    # CHECK REQUIRED COLUMNS
    # =========================
    required = [
        "Total_Abundance",
        "domain",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species"
    ]
    
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
    # =========================
    # BUILD KRONA TABLE
    # =========================
    # Build without dataset first, so the final size is fixed after filtering
    krona = pd.DataFrame(index=df.index)
    
    # amount
    krona["amount"] = pd.to_numeric(df["Total_Abundance"], errors="coerce").fillna(0)
    
    # optional blank column
    if INCLUDE_BLANK_COLUMN:
        krona[""] = ""
    
    # taxonomy columns
    krona["domain"] = df["domain"].fillna("").astype(str)
    krona["phylum"] = df["phylum"].fillna("").astype(str)
    krona["class"] = df["class"].fillna("").astype(str)
    krona["order"] = df["order"].fillna("").astype(str)
    krona["family"] = df["family"].fillna("").astype(str)
    krona["genus"] = df["genus"].fillna("").astype(str)
    krona["species"] = df["species"].fillna("").astype(str)
    
    # =========================
    # REMOVE ZERO-ABUNDANCE ROWS
    # =========================
    before = len(krona)
    krona = krona[krona["amount"] > 0].copy()
    after = len(krona)
    
    print(f"Rows after removing zero-abundance entries: {after:,} (removed {before - after:,})")
    
    # =========================
    # ADD DATASET COLUMN LAST
    # =========================
    # Insert as first column and repeat the dataset name for every row
    krona.insert(0, "dataset", [DATASET_NAME] * len(krona))
    
    # =========================
    # SAVE OUTPUT
    # =========================
    krona.to_csv(OUTPUT_FILE, index=False, encoding="utf-8-sig")
    
    print(f"Saved Krona file: {OUTPUT_FILE.resolve()}")
    
    # =========================
    # PREVIEW
    # =========================
    print("\nPreview of output:")
    print(krona.head())
    
    print("\nColumn names:")
    print(list(krona.columns))

if __name__ == "__main__":
    main()
