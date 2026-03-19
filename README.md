# DB-Builder

> Reproducible metabarcoding database construction pipeline

DB-Builder is a reproducible Snakemake-based pipeline for integrating metabarcoding datasets into standardized biodiversity databases. It combines taxonomic annotations, representative sequences, and abundance tables into unified matrices, and enables taxonomic collapsing of sequence variants using SINTAX confidence thresholds. This facilitates large-scale ecological and biodiversity analyses. Functional annotations for fungi are derived from the [fungaltraits database](https://github.com/traitecoevo/fungaltraits).

---

## 📑 Table of Contents

* [Workflow](#workflow)
* [Installation](#installation)
* [Usage](#usage)
* [Input](#input)
* [Output](#output)
* [Pipeline Details](#pipeline-details)
* [Common Issues](#common-issues)
* [License](#license)

---

## ⚙️ Workflow

The pipeline consists of three main stages followed by accessory stages:

### 🧩 Part 1 — Dataset Integration

* Automatic detection of complete datasets
* Abundance table used as **master reference**
* Integration of:

  * SINTAX taxonomy
  * FASTA sequences
* Generates:

  * Concatenated datasets
  * Fungi subset
  * Non-annotated subset
  * Full eukaryote dataset

---

### 🧬 Part 2 — Taxonomic Collapse & Database Construction

* OTU collapsing based on selected strategy:

  * `species_only`
  * `genus`
  * `all` (deepest confident rank)
* Representative sequence selection (longest sequence with complete taxonomy)
* Abundance aggregation across samples
* Taxonomy sanitization based on confidence thresholds
* SPPN assignment (unique identifiers)
* Integration of fungal ecological traits
* Final database generation
* Generates:

  * Final database

---

### 📦 Part 3 — Export for R (phyloseq-ready)

* Generates:

  * taxonomy table
  * abundance table
  * sequence table
  * fungal traits table

---

### 📦 Part 4 — Build Krona2 file (accessory)

* Generates:

  * Krona2-ready table

---

### 📦 Part 5 — Subset by taxa (accessory/optional)

* Generates:

  * Count table for the specified taxa 
  * Sequences file

---

## 🛠 Installation

### Requirements

* Linux (recommended)
* Python ≥ 3.8
* Snakemake ≥ 7
* Conda (recommended)

### Setup

```bash
git clone https://github.com/Burn121212/DB-Builder.git
cd DB-Builder
```

Create environment:

```bash
conda create -c conda-forge -c bioconda -n dbbuilder snakemake
conda activate dbbuilder
```

(Optional but recommended):

```bash
conda config --set channel_priority strict
```

---

## 🚀 Usage

```bash
snakemake --use-conda --config input_directory="input_directory" output_directory="output_directory"
```

---

### ⚙️ Parameters

| Parameter           | Description                    |
| ------------------- | ------------------------------ |
| `input_directory`   | Path to input datasets         |
| `output_directory`  | DB-Builder working directory   |
| `collapse_mode`     | `fungi` or `all_eukaryotes`    |
| `collapse_strategy` | `species_only`, `genus`, `all` |
| `p-value_threshold` | Confidence threshold (0.0–1.0) |
| `taxonomic_filter`* | Taxon to subset                |
|`taxonmic_rank`*     | Rank of the taxon to subset    |

Default parameters are `all_eukaryotes`, `species_only`, and `1.0`.

*optional

---

### ▶️ Example

Run the workflow:

```bash
snakemake --cores 4 --use-conda --config input_directory="~/my_dataset/dbbuild_input" output_directory="~/my_dataset/dbbuild_out"
```

Run the workflow and retrieve a subset for a specific taxon:

```bash
snakemake --cores 4 --use-conda --config input_directory="~/my_dataset/dbbuild_input" output_directory="~/my_dataset/dbbuild_out" taxonomic_rank="family" taxonomic_filter="Amanitaceae"
```
---

## 📥 Input

Each dataset must include:

```
<PREFIX>_otu_abundance_table.csv
<PREFIX>_sequences.fasta
<PREFIX>_taxonomy.sintax.txt
```

### Requirements

* Matching OTU IDs across all files
* SINTAX taxonomy format
* FASTA headers correspond to OTUs

---

## 📤 Output

### Part 1

```
concatenated_tables/
├── <PREFIX>_concatenated.csv
├── Fungi_concatenated/<PREFIX>_fungi.csv
├── Non_annotated/<PREFIX>_non_annotated.csv
├── AllEuk_concatenated/<PREFIX>_all_eukaryotes.csv
```

---

### Part 2 (Main output)

```
FINAL_DB/
├── Final_Database_<strategy>_all_eukaryotes_p<value>.csv
├── Final_Database_<strategy>_fungi_p<value>.csv
├── Collapse_Log_<strategy>_<suffix>_p<value>.txt
```

---

### Part 3 (For R)

```
For_R/<tag>/
├── taxonomy.csv
├── abundance_table.csv
├── sequences.fasta
├── fungal_traits_table.csv
```

---

### Part 4 (For Krona2)

```
krona2/
├── <dataset_name>_krona_<collapse_strategy>_<collapse_mode>_p<pvalue_threshold>.csv
```

---

### Part 5 

```
subset/
├── <dataset_name>_subset_<taxonomic_rank>_<taxonomic_filter>_p<pvalue_threshold>.csv
├── <dataset_name>_subset_<taxonomic_rank>_<taxonomic_filter>_p<pvalue_threshold>.fasta
```

---

## 🔬 Pipeline Details

### Dataset Discovery

Only complete datasets (abundance + taxonomy + sequences) are processed automatically. 

---

### Master Table Logic

* Abundance tables define the OTU universe
* Only OTUs present in abundance tables are retained

---

### Collapse Strategies

#### `species_only`

* Collapse only if species-level annotation meets threshold

#### `genus`

* Collapse at genus level if confidence threshold is met

#### `all`

* Recursive collapse using deepest confident rank:

```
species → genus → family → order → class → phylum → domain
```

---

### Representative Selection

* Requires **complete taxonomy (d,p,c,o,f,g,s)**
* Selects **longest sequence**

---

### Taxonomy Sanitization

* Removes low-confidence taxonomic assignments
* Applies cascading consistency across ranks

---

### SPPN Assignment

* Unique identifiers based on deepest confident taxon
* Enables consistent downstream analysis

---

### Fungal Traits Integration

* Matches genus names to traits database
* Adds ecological metadata

---

### Subset taxa

* Given a rank (e.g. "family") and a taxon (e.g. "Amanitaceae") or list of taxa separated by a comma (e.g. "Amanitaceae,Boletaceae") provides a table and a sequence file that match the subset parameters `taxonomic_rank` and `taxonomic_filter` 

---

### For_R Export

* Produces phyloseq-compatible files
* Ensures interoperability with R workflows

---

## 📜 License

This project is licensed under:

**Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)**
