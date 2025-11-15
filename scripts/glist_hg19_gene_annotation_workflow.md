## Annotating VCF Variants with Gene Names Using `bcftools` and `glist-hg19`

This guide shows how to add **gene names** to every variant in a VCF by:

1. Using the PLINK **gene range** file `glist-hg19`
2. Preparing it for use with `bcftools annotate`
3. Running `bcftools annotate` to add an `INFO/GENE` field
4. Verifying that gene names were correctly added

The goal is to end up with an annotated VCF where each variant has a `GENE=` tag in the INFO field indicating which gene it falls in.

---

### 1. Requirements and Inputs

#### 1.1 Software you need

You should have the following tools installed and on your `PATH`:

- `bcftools`
- `bgzip` (usually from `htslib` or `tabix` package)
- `tabix`

You can check with:

```bash
bcftools --version
bgzip --help 2>&1 | head -n 1
tabix --version
```

What this does:  
These commands confirm that the tools are installed and available. If any command is “not found”, install the missing tool before continuing.

---

#### 1.2 Files you need

You will need:

- A `VCF file` you want to annotate, compressed and indexed:

```bash
# Example input VCF (change this path to your own)
ls -lh data/cohort.gtc.GRCh37.norm.vcf.gz
ls -lh data/cohort.gtc.GRCh37.norm.vcf.gz.tbi
```

- The PLINK gene range file: `glist-hg19` (we’ll download this below).

---

### 2. Download the PLINK gene range file (glist-hg19)

We will store the gene range file in a reference/gene_ranges_hg19 directory inside your project.

```bash
# Go to your project root (edit this path to match your setup)
cd /path/to/your/project

# Create a directory for gene ranges
mkdir -p reference/gene_ranges_hg19
cd reference/gene_ranges_hg19
```

What this does:

- Creates a dedicated folder for the hg19 gene ranges.
- Moves into that folder so that all following commands run there.

---

#### 2.1 Download glist-hg19 (option A: wget)

If you have the direct URL for glist-hg19, use:

```bash
# Replace the URL with the actual PLINK glist-hg19 link
wget -O glist-hg19 "https://example.org/path/to/glist-hg19"
```

What this does:

- Downloads the gene range file and saves it locally as glist-hg19.

---

#### 2.2 Download glist-hg19 (option B: manual copy)

If you downloaded glist-hg19 in a browser to your local machine, copy it into your project:

```bash
# Still in reference/gene_ranges_hg19
cp /path/to/Downloads/glist-hg19 .
ls -lh glist-hg19
```

What this does:

- Copies glist-hg19 from wherever you saved it (e.g., Downloads) into reference/gene_ranges_hg19.
- Lists the file so you can confirm it’s there.

---

### 3. Inspect and prepare glist-hg19 for bcftools annotate

Our goal is to turn glist-hg19 into a tabix-indexed, bgzipped file that bcftools can query by genomic region.

The format we want is:

- One line per gene
- 4 columns: CHROM START END GENE
- Columns separated by TAB characters (not spaces)
- Sorted by chromosome and start position

---

#### 3.1 Quick inspection

```bash
# Look at the first few lines
head glist-hg19

# Count how many genes/lines
wc -l glist-hg19
```

What this does:

- head shows you the structure of the file (chromosome, start, end, gene name).
- wc -l tells you how many gene records are present.

You will typically see something like:

```sh
1 11873 14409 DDX11L1
1 14362 29370 WASH7P
```

---

#### 3.2 Ensure the file is sorted

We want the file sorted by:

- Column 1: chromosome
- Column 2: start position (numeric)

```bash
# Sort by chromosome and start position
sort -k1,1 -k2,2n glist-hg19 > glist-hg19.sorted

# Replace the original with the sorted version
mv glist-hg19.sorted glist-hg19
```

What this does:

- Creates a new sorted file and then replaces the original with the sorted one, ensuring tabix will work correctly later.

---

#### 3.3 Convert to TAB-delimited format

Some versions of `glist-hg19` may use spaces. tabix requires tab-delimited files. We’ll reconstruct the file with explicit TABs and exactly 4 columns.

```bash
# Rebuild with explicit TAB separators between columns
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' glist-hg19 > glist-hg19.tab

# Replace original file with the tab-delimited version
mv glist-hg19.tab glist-hg19
```

What this does:

- Uses awk to print the first 4 columns with TABs between them.
- Ensures a clean, 4-column, TAB-separated file for indexing.

You can confirm TABs are present (shown as ^I) with:

```bash
head glist-hg19 | cat -t
```

You should now see lines like:

`1^I11873^I14409^IDDX11L1`

---

#### 4. Compress and index glist-hg19 for region-based lookup

bcftools annotate expects the annotation file to be bgzipped and tabix-indexed so it can quickly find which gene overlaps each variant.

---

#### 4.1 Compress glist-hg19 with bgzip

```bash
# Compress; -f overwrites any existing .gz file
bgzip -f glist-hg19

# List files to confirm
ls -lh
```

You should see a glist-hg19.gz file.

What this does:

- bgzip compresses the file in a way that is compatible with tabix and bcftools.

(Optional) Check inside the compressed file:

```bash
zcat glist-hg19.gz | head
zcat glist-hg19.gz | wc -l
```

---

### 4.2 Index glist-hg19.gz with tabix

```bash
# -s 1: chromosome is column 1
# -b 2: start position is column 2
# -e 3: end position is column 3
tabix -s 1 -b 2 -e 3 glist-hg19.gz

# Confirm both files exist
ls
```

you should now see:

- `glist-hg19.gz`
- `glist-hg19.gz.tbi`

What this does:

- Creates an index file (.tbi) that allows fast random access to gene intervals, enabling bcftools to find which gene(s) overlap a given variant position.

If you get an error (e.g., “Failed to parse TBX_GENERIC”), double-check that the file is TAB-delimited and has at least 3 columns; then repeat Step 3.3 and Step 4.

---

### 5. Create a VCF header snippet for the GENE INFO field

VCF files must define each INFO tag in the header. We will create a small file that describes the new GENE field.

```bash
# Stay in reference/gene_ranges_hg19
pwd

# Create gene_header.txt with a single INFO definition
cat << 'EOF' > gene_header.txt
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name from PLINK glist-hg19 annotation">
EOF

# Inspect the file
cat gene_header.txt
```

What this does:

- Writes a one-line header definition for an INFO field called GENE.
- Number=1 means there is one gene name per variant (for variants overlapping multiple genes, you may need to adapt this logic later).

---

### 6. Annotate your VCF with gene names using bcftools annotate

Now we have:

- `glist-hg19.gz` and `glist-hg19.gz.tbi` (annotation source)
- `gene_header.txt` (INFO definition)

We can use `bcftools annotate` to add a `GENE=` tag to each variant.

Assume your input VCF is:
`data/cohort.gtc.GRCh37.norm.vcf.gz`

---

#### 6.1 Run bcftools annotate

```bash
cd /path/to/your/project

bcftools annotate \
  -a reference/gene_ranges_hg19/glist-hg19.gz \
  -h reference/gene_ranges_hg19/gene_header.txt \
  -c CHROM,FROM,TO,GENE \
  -Oz -o data/cohort.gtc.GRCh37.annotated.vcf.gz \
  data/cohort.gtc.GRCh37.norm.vcf.gz
```

**What each option does:**

- -a reference/gene_ranges_hg19/glist-hg19.gz  
   Tells bcftools which annotation file to use (our gene ranges).
- -h reference/gene_ranges_hg19/gene_header.txt  
   Adds the ##INFO=<ID=GENE,...> line into the output VCF header.
- -c CHROM,FROM,TO,GENE  
   Maps columns in the annotation file to VCF fields:

   - CHROM → chromosome
   - FROM → start coordinate
   - TO → end coordinate
   - GENE → the value for the INFO/GENE field
      Internally, bcftools uses the region (CHROM:FROM-TO) to decide if a variant falls in that gene.

- -Oz
   Writes the output as a bgzipped VCF (.vcf.gz).
- -o data/cohort.gtc.GRCh37.annotated.vcf.gz
   Path for the annotated output VCF.
- Final argument data/cohort.gtc.GRCh37.norm.vcf.gz
   Your original, unannotated VCF.

---

#### 6.2 Index the annotated VCF

After annotation, index the new VCF so that tools can query it efficiently:

```bash
tabix -p vcf data/cohort.gtc.GRCh37.annotated.vcf.gz
```

What this does:

- -p vcf tells tabix that this is a VCF file.
- Creates data/cohort.gtc.GRCh37.annotated.vcf.gz.tbi so the annotated VCF can be used in downstream tools.

---

### 7. Verify that gene annotation worked

Here we perform a few simple checks to confirm that:

- The GENE field is present in the header.
- Variants have GENE= entries in their INFO field.
- We can export a simple table: CHROM POS ID GENE REF ALT.

---

#### 7.1 Check the header for the GENE INFO definition

```bash
bcftools view -h data/cohort.gtc.GRCh37.annotated.vcf.gz | grep "INFO=<ID=GENE"
```

What this does:

- Shows the header line that defines the GENE field.
- Confirms that our gene_header.txt was successfully added to the VCF.

You should see something like:

```sh
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name from PLINK glist-hg19 annotation">
```

#### 7.2 Look at a few annotated variants

```bash
bcftools view data/cohort.gtc.GRCh37.annotated.vcf.gz | head -n 20
```

What this does:

- Prints the first 20 lines of the VCF (excluding the header).
- Look in the INFO field for GENE=... next to each variant.

An example INFO field might look like:
`INFO=AC=1;AF=0.003;AN=312;GENE=DDX11L1`

---

### 7.3 Export a simple variant–gene table (TSV or CSV)

You can also convert the annotated VCF into a text table for easier inspection or downstream analysis.

**TSV (tab-separated):**

```bash
bcftools query \
  -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%GENE\n' \
  data/cohort.gtc.GRCh37.annotated.vcf.gz \
  > data/cohort.gtc.GRCh37.annotated.gene_table.tsv
```

**CSV (comma-separated):**

```bash
bcftools query \
  -f '%CHROM,%POS,%ID,%REF,%ALT,%GENE\n' \
  data/cohort.gtc.GRCh37.annotated.vcf.gz \
  > data/cohort.gtc.GRCh37.annotated.gene_table.csv
```

What this does:

- Uses bcftools query to extract specific fields from each variant:
   - Chromosome (%CHROM)
   - Position (%POS)
   - Variant ID (%ID)
   - Reference allele (%REF)
   - Alternate allele (%ALT)
   - Gene name from the annotation (%GENE)

- Saves them into a plain text table that you can open in R, Python, Excel, etc.

You can quickly preview the table with:

```bash
head data/cohort.gtc.GRCh37.annotated.gene_table.tsv
```

---

### 8. Summary

By following this workflow, you have:

1. Downloaded and cleaned the PLINK glist-hg19 gene range file.
2. Converted it to a TAB-delimited, bgzipped, and tabix-indexed annotation source.
3. Created a VCF header snippet defining a GENE INFO field.
4. Used bcftools annotate to attach gene names to each variant in your VCF.
5. Verified the annotation and exported a variant–gene table.

You can now use the annotated VCF and the exported tables for downstream analyses such as:

- Filtering variants by gene
- Linking variants to gene-level functional data
- Integrating with expression or pathway analysis

---

### Optional: One-shot script to prepare `glist-hg19` and header

If you prefer not to run each preparation step manually, this repository includes a helper script that performs all steps automatically:

```bash
bash scripts/glist_hg19_gene_annotation_prepare.sh

This script performs exactly the same operations described in this document, but in a fully automated way:
- Creates the directory `reference/gene_ranges_hg19` and switches into it  
    This corresponds to the “make directory and move into it” step.

- Downloads `glist-hg19` from the official PLINK resources (if not already present)  
    Matches the “download via wget or manual copy” instructions.

- Checks whether `glist-hg19` is sorted and sorts by chromosome and start position if required  
    Mirrors the “ensure the file is sorted” step.

- Rebuilds the file as a clean, TAB-delimited 4-column format (`CHR`, `START`, `END`, `GENE`)  
    Same as the “convert to TAB-delimited format” step.

- Compresses the file with `bgzip` and builds a `tabix` index (`-s 1 -b 2 -e 3`)  
    Exactly as described in the “compress and index” section.

- Creates a `gene_header.txt` file defining the VCF `INFO/GENE` annotation definition  
    Same as the “create a VCF header snippet” instructions.

- Prints an example `bcftools annotate` command for convenience  
    Matches the annotation example provided at the end of this guide.

After running this script, you will have the following ready-to-use files:
- `reference/gene_ranges_hg19/glist-hg19.gz`
- `reference/gene_ranges_hg19/glist-hg19.gz.tbi`
- `reference/gene_ranges_hg19/gene_header.txt`
```