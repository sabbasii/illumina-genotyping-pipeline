## Exporting SNP–Gene–Genotype Tables from an Annotated VCF

> ### 🔗 Related Document
>
> To prepare `glist-hg19` and run the gene-annotation step, see:  
> __`glist_hg19_gene_annotation_workflow.md`__

This document shows how to:

1. Confirm that your VCF already contains **gene names** in `INFO/GENE`
2. (If needed) **index** the annotated VCF
3. Export a **TSV** and **CSV** table with columns:  
   `ID | CHROM | POS | REF | ALT | GENE | Genotypes…`

This guide assumes the **gene annotation step is already done** (using `bcftools annotate` and `glist-hg19`) as described in:

- `glist_hg19_gene_annotation_workflow.md`

---

### 0. Assumptions and Required Inputs

We assume you already have:

- An **annotated VCF** (GENE field present):

```text
output/genotype_run1/vcf/cohort.gtc.GRCh37.annotated.vcf.gz

```

The corresponding index file (if not, we’ll create it):

```text
output/genotype_run1/vcf/cohort.gtc.GRCh37.annotated.vcf.gz.tbi
```

The annotated VCF should already contain an INFO/GENE field produced by the previous workflow; we won’t repeat that annotation step here.

---

### 1. Check and (If Needed) Create the Index for the Annotated VCF

Move into the VCF directory:

```bash
cd ~/git_projects/illumina-genotyping-pipeline/output/genotype_run1/vcf
```

Check whether the index exists:

```bash
ls cohort.gtc.GRCh37.annotated.vcf.gz.tbi
```

- If the file exists, you can skip indexing.
- If it does not exist, create it:

```bash
tabix -p vcf cohort.gtc.GRCh37.annotated.vcf.gz
```

What this does:

- `tabix -p vcf` builds a `.tbi` index for the annotated VCF, which is required for efficient region-based queries and some downstream tools.

---

### 2. Verify That the GENE Annotation Is Present

Before exporting tables, quickly confirm that the VCF contains an INFO/GENE field.

#### 2.1 Check the header for GENE

```bash
bcftools view -h cohort.gtc.GRCh37.annotated.vcf.gz | grep 'INFO=<ID=GENE'
```

You should see a line like:

```text
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name from PLINK glist-hg19 annotation">
```

What this does:

- Confirms that the `GENE` INFO field is defined in the header, meaning the previous annotation step succeeded.

---

#### 2.2 Spot-check a few variants with GENE values

```basg
bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/GENE\n' \
  cohort.gtc.GRCh37.annotated.vcf.gz | head
```

Example output (simplified):

```text
rs9701055    1    565433    A    G    .
rs9651229    1    567667    C    T    .
GSA-rs200599638 1 752918    G    A    FAM87B
rs12127425   1    794332    T    C    LINC01128
...
```

- `.` in the last column → the variant does not overlap any gene interval (in this list).
- A gene symbol like FAM87B or LINC01128 → successfully annotated.

---

### 3. Export a SNP–Gene–Genotype Table as TSV

We will create a wide table with:

```text
ID  CHROM  POS  REF  ALT  GENE  SAMPLE1_GT  SAMPLE2_GT  ...
```

Use `bcftools query` to output a tab-separated file:

```bash
bcftools query -H \
  -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/GENE[\t%GT]\n' \
  cohort.gtc.GRCh37.annotated.vcf.gz \
  > cohort_snps_genes_genotypes.tsv
```

#### 3.1 Format explanation

- `-H`
   Adds a header row with column names.
- `-f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/GENE[\t%GT]\n'`
   Tells bcftools which fields to print:

   - %ID → variant ID (e.g., rsID or array ID)
   - %CHROM → chromosome
   - %POS → genomic position
   - %REF, %ALT → reference and alternate alleles
   - %INFO/GENE → gene name from the annotation step
   - [\t%GT] → for each sample in the VCF, add a TAB and that sample’s genotype (e.g., 0/0, 0/1, 1/1)

- `> cohort_snps_genes_genotypes.tsv`
   Saves the output as a tab-separated values (TSV) file.

Preview the first few lines:

```bash
head cohort_snps_genes_genotypes.tsv | sed -n '1,5p'
```

---

### 4. Export a SNP–Gene–Genotype Table as CSV

If you prefer a comma-separated file (CSV) for Excel or other tools, you can either:

1. Export directly as CSV with `bcftools`, or
2. Convert the TSV to CSV.

#### 4.1 Direct CSV export via `bcftools query`

```bash
bcftools query -H \
  -f '%ID,%CHROM,%POS,%REF,%ALT,%INFO/GENE[%GT,]\n' \
  cohort.gtc.GRCh37.annotated.vcf.gz \
  > cohort_snps_genes_genotypes.csv
```

What this does:

- Uses commas instead of tabs.
- `[%GT,]` prints each sample’s genotype followed by a comma.
- Produces a wide CSV where columns after `GENE` correspond to sample genotypes.

---

#### 4.2 Convert an existing TSV to CSV

If you already created cohort_snps_genes_genotypes.tsv, you can convert it:

```bash
tr '\t' ',' < cohort_snps_genes_genotypes.tsv \
  > cohort_snps_genes_genotypes.csv
```

What this does:

- Reads the TSV file.
- Replaces each TAB (`\t`) with a comma.
- Writes the result as a CSV file.

---

### 5. Summary of Outputs

After running the commands in this workflow, you should have:

- `output/genotype_run1/vcf/cohort.gtc.GRCh37.annotated.vcf.gz`
   → VCF with `INFO/GENE` for each variant (produced in the previous workflow).
- `output/genotype_run1/vcf/cohort.gtc.GRCh37.annotated.vcf.gz.tbi`
   → Tabix index for the annotated VCF.
- `output/genotype_run1/vcf/cohort_snps_genes_genotypes.tsv`
   → Wide TSV table: `ID, CHROM, POS, REF, ALT, GENE, genotypes` for all samples.
- `output/genotype_run1/vcf/cohort_snps_genes_genotypes.csv`
   → Same table in CSV format for Excel / R / Python.

These tables are ready for downstream analysis, such as:

- Filtering variants by gene
- Per-gene or per-region summaries
- Integrating with expression, pathway, or clinical data

---

### Optional: One-shot script to generate SNP–Gene–Genotype tables

If you prefer not to run each command manually, this repository includes a helper script that runs the full export workflow in one go:

```bash
bash scripts/vcf_gene_annotation_export.sh
```

The script:

- Ensures you are working in output/genotype_run1/vcf
- Checks that cohort.gtc.GRCh37.annotated.vcf.gz exists
- Creates a tabix index (cohort.gtc.GRCh37.annotated.vcf.gz.tbi) if it is missing
- Verifies that the VCF header contains an INFO/GENE definition
- Performs a small spot-check of variants to show the GENE field
- Exports a wide SNP–gene–genotype table as:

   - cohort_snps_genes_genotypes.tsv
   - cohort_snps_genes_genotypes.csv (converted from the TSV)

The outputs of vcf_gene_annotation_export.sh match the tables described in this document and are ready for downstream analysis.

---

### Optional: R script for quick CSV inspection

For an interactive QC / exploration of the exported CSV, you can use:

```bash
Rscript scripts/inspect_annotated_snps.R

This script:

- Reads cohort_snps_genes_genotypes.csv
- Reports basic structure (rows, columns, column names, types)
- Summarizes how many SNPs have gene annotations and how many genes are represented
- Shows genotype distributions for a few samples
- Provides an example of extracting all SNPs for a given gene (default: FAM87B)