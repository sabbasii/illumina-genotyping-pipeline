# VCF/BCF Quality Checks

This page explains the routine checks we run after converting Illumina GTC files into VCF/BCF format with `bcftools +gtc2vcf`. These QC steps help confirm that the conversion succeeded, the file is valid, and the data look reasonable before downstream analysis.

---

## Files
- **Input**:  
  `output/<RUN>/vcf/cohort.gtc.<REF_BUILD>.norm.vcf.gz`  
  (or the corresponding `.bcf` file if you prefer binary)
- **Outputs**:  
  `output/<RUN>/qc/*`  
  (small text files & plots safe to version-control)

---

## Checks & Why

### 1. Header & contigs  
- **Files**: `header.txt`, `contigs.list`  
- **What it does**: Dumps the VCF header and extracts all `##contig` lines.  
- **Why**: The header contains metadata (reference genome build, command line, bcftools version). Contig names must match your reference FASTA (e.g., `1..22, X, Y` for GRCh37 vs `chr1..chr22, chrX, chrY` for GRCh38). Mismatches cause downstream errors in imputation, annotation, or liftover.  
- **Interpretation**: Contigs listed here should exactly match your `$REFERENCE_FASTA`. If not, fix before downstream steps.

**Note:** If the VCF header lacks `##contig=` lines, `qc_vcf.sh` automatically injects them from the FASTA index (`reference.fa.fai`) via `bcftools reheader`. This is a header-only change; it does not modify genotypes or records.

---

### 2. Sample count  
- **Files**: `samples.list`, `samples.count`  
- **What it does**: Lists all sample IDs in the VCF header.  
- **Why**: Confirms that every `.gtc` processed shows up in the VCF. The number must equal the count of `.gtc` files detected in the run.  
- **Interpretation**: If you expect 288 samples and only see 287, check logs for failed samples.

---

### 3. Variant count  
- **File**: `variants.count`  
- **What it does**: Uses the VCF index to count total variant records.  
- **Why**: Sanity check against the manifest (~700–750k SNPs for GSA v3). Large deviations may indicate missing or duplicate records.  
- **Interpretation**: Count should be close to manifest expectations. A big mismatch often means wrong manifest or reference.

---

### 4. bcftools stats + plots  
- **Files**: `bcftools.stats.txt`, `vcfstats_plots/`  
- **What it does**: Generates summary metrics (Ti/Tv ratio, singleton rate, multiallelic proportion, depth distribution, etc.) and optional plots.  
- **Why**: Provides a quick snapshot of data quality. Many pipelines expect Ti/Tv ≈ 2.0–2.1 in human data; much lower can indicate contamination or reference mismatch.  
- **Interpretation**:  
  - **Ti/Tv**: ~2 is normal; <1.5 may be suspicious.  
  - **Singleton rate**: High rate can suggest genotyping errors.  
  - **Depth (DP)**: Should cluster around expected array depth (not sequencing, but probe-based).  
  - **Multiallelics**: Should be present but not excessive.

---

### 5. Variants by contig  
- **File**: `variants_by_contig.txt`  
- **What it does**: Counts variants per chromosome.  
- **Why**: Ensures coverage matches array design. If one contig has zero variants, it may not be included in the manifest.  
- **Interpretation**: All standard chromosomes should appear with expected counts; spurious contigs (e.g., `GL000xxx`) might indicate manifest artifacts.

---

### 6. REF/ALT sanity checks  
- **Files**: `refcheck.log`, `ref_mismatch.summary`  
- **What it does**: Uses `bcftools norm -c ws` to check for REF/ALT mismatches against the reference FASTA.  
- **Why**: If the REF allele in the VCF doesn’t match the reference genome, downstream tools may fail or silently misinterpret variants.  
- **Interpretation**: Ideally, the summary is empty. If mismatches are found, verify manifest genome build or reference FASTA.

---

### 7. Allele frequency snapshot  
- **File**: `af.tsv.gz`  
- **What it does**: Extracts allele frequencies (AF) from INFO tags (populated by `bcftools +fill-tags`).  
- **Why**: Gives a quick distribution of allele frequencies across the cohort. Useful for sanity checks (e.g., very rare variants vs common SNPs).  
- **Interpretation**: Expect a mix of rare and common variants. Uniform or extreme distributions can indicate issues with calling.

**### Note on AF availability  
Your normalized VCF may not contain `INFO/AF` by default. Our QC script computes AF on the fly using:  bcftools +fill-tags -t AF | bcftools query ...

- This injects AF (and the proper header line) into the stream without rewriting your original file.  
- If you prefer to persist AF in a final deliverable, you can create `cohort.gtc.<REF_BUILD>.norm.withAF.vcf.gz` by piping `+fill-tags` into `bcftools view -Oz` and indexing it.**  
---

### 8. Primary-only deliverable (using `MAKE_PRIMARY=1`)
RUN:
```bash
MAKE_PRIMARY=1 bash scripts/qc_vcf.sh
```
to generate a clean, primary-only VCF deliverable for downstream analyses, while keeping your full normalized VCF as the source of truth.

#### What this is
`MAKE_PRIMARY` is an optional toggle for `scripts/qc_vcf.sh`. When enabled, the script writes an **extra** VCF that contains **only the primary chromosomes**.  
- **Primary chromosomes (GRCh37/hg19):** `1..22, X, Y, MT`  
- **Primary chromosomes (GRCh38 with chr prefixes):** `chr1..chr22, chrX, chrY, chrM`  

Your **original normalized VCF** remains unchanged.

### Why use it
- **Cleaner downstream pipelines:** Many GWAS/QC/imputation workflows ignore unplaced scaffolds (`GL0002xx`, etc.). A primary-only VCF avoids surprises and mismatches.  
- **Better tool compatibility:** Some tools/servers expect only primary contigs by default.  
- **Smaller, faster files:** Fewer contigs → faster operations, simpler sharing.  
- **Reproducibility:** Produces a standard, “ready-for-GWAS” deliverable while keeping the full VCF intact.


#### How to run it
From the project root:
```bash
source scripts/00_config.sh
MAKE_PRIMARY=1 bash scripts/qc_vcf.sh
```
### What you’ll get
A new, bgzipped & indexed VCF alongside your normalized file:
    output/<RUN>/vcf/cohort.gtc.<REF_BUILD>.norm.primary.vcf.gz
    output/<RUN>/vcf/cohort.gtc.<REF_BUILD>.norm.primary.vcf.gz.tbi

### What the script does under the hood
1. Builds the primary contig list for your reference (e.g., 1..22, X, Y, MT for GRCh37).
2. Joins that list into a bcftools regions string (e.g., 1,2,...,22,X,Y,MT).
3. Subsets your normalized VCF:
    bcftools view -r 1,2,...,22,X,Y,MT -Oz -o *.primary.vcf.gz
    bcftools index -t *.primary.vcf.gz
4. Leaves the original VCF untouched.

### When to use (and when not to)
- Use it when sharing data for GWAS/QC/imputation or feeding pipelines that don’t want scaffolds/alternates.
- Skip it if your analysis explicitly needs variants on non-primary contigs (e.g., unplaced scaffolds).

---

### 9. PLINK-based QC  
- **Files**: `plink/cohort.imiss`, `cohort.lmiss`, `cohort.afreq`, `cohort.hardy`  
- **What it does**: Runs PLINK2 QC modules for missingness, allele frequency, and Hardy–Weinberg equilibrium.  
- **Why**: Provides per-sample and per-variant QC metrics widely used in GWAS pipelines.  
- **Interpretation**:  
  - **.imiss**: High missingness samples may be outliers.  
  - **.lmiss**: Variants with high missingness should be excluded.  
  - **.hardy**: Variants with extreme HWE deviation can be problematic.  
  - **.afreq**: Allele frequency spectrum should look plausible (compare with 1000 Genomes if desired).

---

## PLINK QC (simple guide)

PLINK is a tool that checks if our genotype data looks healthy. After we create the normalized VCF, the QC script runs PLINK2 to compute:

- **Missingness**: how much data is missing per **sample** (`.smiss.gz`) and per **variant** (`.vmiss.gz`)
- **Allele frequencies**: how common each allele is (`.afreq.gz`)
- **Hardy–Weinberg equilibrium (HWE)**: population check for each variant (`.hardy.gz`)

### What the script runs
- **Autosomes (1–22)** always → avoids chrX issues when sex is unknown.
- **chrX (optional)** → only if a PSAM with SEX is available (see below).

### Making chrX possible (linking sex from the sample sheet)
We generate a **PSAM** file from the Illumina sample sheet and feed it to PLINK:

1. Build PSAM with sex:
```bash
   bash scripts/make_psam_from_sample_sheet.sh \
     input_data/sample_sheet/<your_sample_sheet>.csv \
     metadata/cohort.sex.psam
```
2. The QC script picks it up automatically via PSAM_SEX (set in scripts/00_config.sh).

With this in place, PLINK runs chrX with --psam and the proper pseudoautosomal region split (--split-par hg19 for GRCh37, hg38 for GRCh38).

### Where the files go and what they mean
Location: output/<RUN>/qc/plink/

    - **cohort.smiss.gz →** Sample missingness. High missingness samples may be excluded.
    - **cohort.vmiss.gz →** Variant missingness. Sites with lots of missing data may be dropped.
    - **cohort.afreq.gz →** Allele frequency. Helps spot odd distributions and compute MAF.
    - **cohort.hardy.gz →** Hardy–Weinberg. Strong deviations can indicate genotyping problems.
    - **cohort.chrX.*.gz →** chrX-only versions of the same stats (present only when PSAM with sex is provided).
    - **cohort.log →** What PLINK did (helpful for troubleshooting).

Quick checks
```bash
    zcat output/<RUN>/qc/plink/cohort.smiss.gz | head
    zcat output/<RUN>/qc/plink/cohort.vmiss.gz | head
    zcat output/<RUN>/qc/plink/cohort.afreq.gz | head
    zcat output/<RUN>/qc/plink/cohort.hardy.gz | head
```
Tip: Keep your filtering thresholds simple and visible in your analysis notebooks or a FILTERS.md (e.g., sample F_MISS <= 0.02, site F_MISS <= 0.02–0.05, HWE p >= 1e-6 in controls).

---

## Documentation

- [Setup Instructions](docs/SETUP.md)
- [Output Files](docs/OUTPUTS.md)
- [QC Checks](docs/QC.md)