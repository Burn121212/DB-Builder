# =========================================================
# MODULES
# =========================================================

import subprocess
import glob
import os
import pathlib
from collections import defaultdict

# =========================================================
# FUNCTIONS
# =========================================================

def discover_datasets(input_dir: Path) -> dict:
    """Detect complete datasets in input/ (abundance + sequences + taxonomy)."""
    datasets = defaultdict(dict)

    for file_path in input_dir.glob("*"):
        filename = file_path.name

        if filename.endswith("_otu_abundance_table.csv"):
            prefix = filename.replace("_otu_abundance_table.csv", "")
            datasets[prefix]["abundance"] = file_path

        elif filename.endswith("_sequences.fasta"):
            prefix = filename.replace("_sequences.fasta", "")
            datasets[prefix]["sequences"] = file_path

        elif filename.endswith("_taxonomy.sintax.txt"):
            prefix = filename.replace("_taxonomy.sintax.txt", "")
            datasets[prefix]["taxonomy"] = file_path

    complete = {}
    required = {"abundance", "sequences", "taxonomy"}
    for prefix, files in datasets.items():
        missing = required - set(files.keys())
        if not missing:
            complete[prefix] = files
            print(f"✅ Complete dataset detected: {prefix}")
        else:
            print(f"⚠️ Incomplete dataset: {prefix} - missing: {missing}")
    return list(complete.keys())

# =========================================================
# PARSE COMMAND AND SET VARIABLES
# =========================================================

# minimal required arguments
querydir = Path(config["input_directory"]).expanduser().resolve()
outdir = Path(config["output_directory"]).expanduser().resolve()
LOCBASE = discover_datasets(querydir)

# optional arguments with default options
# for collapse
collapse_mode = config.get('collapse_mode', 'all_eukaryotes')
collapse_strategy = config.get('collapse_strategy', 'species_only')
pvalue_threshold = float(config.get('p-value_threshold', 1.0))

# optional arguments with default options
# for subset by taxonomy
taxonomic_rank = config.get('taxonomic_rank', 'species')
taxonomic_filter = config.get('taxonomic_rank', 'None') # cannot be none
remove_prefix = config.get('remove_prefix', 'None')
tax_filter_mode = config.get('tax_filter_mode', 'exact')
dataset_name = config.get('dataset_name', 'concatenated_DB')

# optional arguments with default options
# for build krona
include_blank_column = config.get('include_blank_column', '1')

# =========================================================
# RULE ALL
# =========================================================

rule all:
    input:
        str(outdir) + f"/FINAL_DB/Final_Database_{collapse_strategy}_{collapse_mode}_p{pvalue_threshold}.csv", # collapse
        str(outdir) + f"/FINAL_DB/Collapse_Log_{collapse_strategy}_{collapse_mode}_p{pvalue_threshold}.txt", # collapse
        expand(str(outdir) + "/concatenated_tables/{lbase}_concatenated.csv", lbase=LOCBASE), # collapse
        str(outdir) + f"/krona2/{dataset_name}_krona_{collapse_strategy}_{collapse_mode}_p{pvalue_threshold}.csv" # krona

# =========================================================
# RULES
# =========================================================

rule run_collapse:
    """ dbbuilder_concatenate.py README:
ATLASMX ITS — DB BUILDER V9.1

Usage:
dbbuilder_concatenation.py <input_directory> <output_prefix> <collapse_mode> <collapse_strategy> <p-value_threshold>

Input directory:
    string indicating the path of the input data, for example "~/atlas/dbbuilder/to_build/input"
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
python snakes/dbbuilder_concatenate.py ~/atlas/dbbuilder/to_build/input/ ~/atlas/dbbuilder/db_test all_eukaryotes species_only 0.8

Details:
PART 1 (Concatenation):
1) Abundance table is MASTER (only OTUs present there are kept).
2) Join taxonomy (SINTAX) + sequences (FASTA).
3) Save:
   - concatenated_tables/<PREFIX>_concatenated.csv
   - concatenated_tables/Fungi_concatenated/<PREFIX>_fungi.csv
   - concatenated_tables/Non_annotated/<PREFIX>_non_annotated.csv
   - concatenated_tables/AllEuk_concatenated/<PREFIX>_all_eukaryotes.csv   (UNFILTERED full table)

PART 2 (Taxonomic collapse + final DB):
4) Collapse according to COLLAPSE_STRATEGY:
   - "species_only" -> collapse only when species p >= P_VALUE_THRESHOLD
   - "genus"        -> collapse only when genus p >= P_VALUE_THRESHOLD
   - "all"          -> recursive lowest-rank collapse using the deepest available rank
                       with p >= P_VALUE_THRESHOLD among:
                       species -> genus -> family -> order -> class -> phylum -> domain
5) Representative is chosen ONLY among rows with COMPLETE taxonomy (d,p,c,o,f,g,s);
   among those, the longest sequence is kept.
6) Sum abundances across sites.
7) Compute Total_Abundance, sort, and build binomial species names.
8) SANITIZE taxonomy columns by p >= SPPN_P_THRESHOLD
   (fix SINTAX cases where names remain despite low confidence).
9) Add fungal lifestyles (Fungal_Traits_DB) after sanitization.
10) Assign SPPN (deepest rank with p >= SPPN_P_THRESHOLD), prune zeros.
11) Save final DB CSVs with strategy and p-value in the filename:
    - FINAL_DB/Final_Database_<strategy>_all_eukaryotes_p{P_VALUE_THRESHOLD}.csv   (when MODE=all_eukaryotes)
    - FINAL_DB/Final_Database_<strategy>_fungi_p{P_VALUE_THRESHOLD}.csv             (subset from all-euk OR main when MODE=fungi)

PART 3 (For_R export):
12) Create For_R/ folder with phyloseq-ready files for each final DB produced:
    For_R/<tag>/taxonomy.csv
    For_R/<tag>/abundance_table.csv
    For_R/<tag>/sequences.fasta
    For_R/<tag>/fungal_traits_table.csv

NOTE:
- MODE controls which tables are used for collapse in Part 2.
- Regardless of MODE, if MODE="all_eukaryotes" we also export a fungi-only final DB + its For_R package.
- This version preserves the original conservative logic for representative selection:
  a group is collapsed only if at least one row in that group has COMPLETE taxonomy.
  Otherwise, the group is left uncollapsed.
- For_R folder names are strategy-specific and threshold-specific to avoid overwriting.
    """
    conda:
        "snakes/pandas.yaml"
    input:
        expand(str(querydir) + "/{lbase}_otu_abundance_table.csv", lbase=LOCBASE),
        expand(str(querydir) + "/{lbase}_sequences.fasta", lbase=LOCBASE),
        expand(str(querydir) + "/{lbase}_taxonomy.sintax.txt", lbase=LOCBASE)
    output:
        str(outdir) + f"/FINAL_DB/Final_Database_{collapse_strategy}_{collapse_mode}_p{pvalue_threshold}.csv", # collapse
        str(outdir) + f"/FINAL_DB/Collapse_Log_{collapse_strategy}_{collapse_mode}_p{pvalue_threshold}.txt", # collapse
        expand(str(outdir) + "/concatenated_tables/{lbase}_concatenated.csv", lbase=LOCBASE) # collapse
    params:
        ' '.join([str(querydir),
        str(outdir),
        collapse_mode,
        collapse_strategy,
        f'{pvalue_threshold}'])
    shell:
        """
        python snakes/dbbuilder_concatenate.py {params}
        """


rule run_krona:
    """ dbbuilder_build_krona.py README
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
python snakes/dbbuilder_build_krona.py ~/atlas/dbbuilder/db_test FINAL_DB/Final_Database_species_only_fungi_p0.8.csv Krona_format_species_only_fungi_p0.8_zOTU.csv fungi_zOTU_p0.8 1
    """
    conda:
        "snakes/pandas.yaml"
    input:
        str(outdir) + f"/FINAL_DB/Final_Database_{collapse_strategy}_{collapse_mode}_p{pvalue_threshold}.csv" # collapse
    output:
        str(outdir) + f"/krona2/{dataset_name}_krona_{collapse_strategy}_{collapse_mode}_p{pvalue_threshold}.csv" # krona
    params:
        ' '.join([str(outdir),
        f"FINAL_DB/Final_Database_{collapse_strategy}_{collapse_mode}_p{pvalue_threshold}.csv",
        f"{dataset_name}_krona_{collapse_strategy}_{collapse_mode}_p{pvalue_threshold}.csv",
        dataset_name,
        include_blank_column])
    shell:
        """
        python snakes/dbbuilder_build_krona.py {params}
        """


rule run_subset:
    """ dbbuilder_subset_by_taxonomy.py README
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

Details:
INPUT columns expected (example):
SPPN, OTU_XX, Original_OTUID, domain..species, *_pvalue,
primary_lifestyle, Secondary_lifestyle, sequence, sequence_lenght, Total_Abundance,
Sitio_AtlasEXP01, Sitio_AtlasEXP02, ...

OUTPUT columns desired:
[fixed columns WITHOUT input lifestyles] + [renamed sites] + [FULL fungal traits columns]
where primary_lifestyle & Secondary_lifestyle come ONLY from traits.
    """
    conda:
        "snakes/pandas.yaml"
    input:
        str(outdir) + f"/FINAL_DB/Final_Database_{collapse_strategy}_{collapse_mode}_p{pvalue_threshold}.csv" # collapse
    output:
        str(outdir) + f"/subset/{dataset_name}_subset_{taxonomic_rank}_{taxonomic_filter}_p{pvalue_threshold}.csv", # subset 
        str(outdir) + f"/subset/{dataset_name}_subset_{taxonomic_rank}_{taxonomic_filter}_p{pvalue_threshold}.fasta" # subset 
    params:
        ' '.join([str(outdir),
        f"FINAL_DB/Final_Database_{collapse_strategy}_{collapse_mode}_p{pvalue_threshold}.csv",
        taxonomic_rank,
        taxonomic_filter,
        f'{dataset_name}_subset_{taxonomic_rank}_{taxonomic_filter}_p{pvalue_threshold}',
        f'{dataset_name}_subset_{taxonomic_rank}_{taxonomic_filter}_p{pvalue_threshold}',
        remove_prefix,
        tax_filter_mode])
    shell:
        """
        python snakes/dbbuilder_subset_by_taxonomy.py {params}
        """
