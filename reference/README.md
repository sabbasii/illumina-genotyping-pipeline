# Reference Genome Guide

This document explains **which reference genome to use**, **where it lives in this repo**, and **how to validate it** for the Illumina genotyping-array pipeline.

---

## 1. Where the reference genome lives

**Note:**  
Reference genomes live in this `reference/` directory and are **ignored by git**.

This is intentional:
- Reference FASTA files are large
- Different users may require different genome builds
- GitHub should not track binary reference files

Each user is responsible for downloading the correct reference locally.

---

## 2. Which reference genome should I use?

### Always check the Illumina manifest

The **single source of truth** is the **`GenomeBuild`** column in the **manifest file you are actually using**.

For this project:
- The manifest indicates **GenomeBuild = 37**
- Therefore, the pipeline must use a **GRCh37-compatible reference**

> Do not assume GRCh38 unless the manifest explicitly says so.

---

## 3. GSAMD and custom array context

- **GSAMD** indicates a **custom Illumina genotyping array**
- Custom arrays are still tied to a **specific reference genome**
- In this project, GSAMD arrays are built on **GRCh37**

This means:
- Probe coordinates
- Allele definitions
- Variant positions  

are all defined relative to **GRCh37**.

---

## 4. Important note on custom arrays (very important)

If your array is **custom**, the **standard Infinium Global Screening Array v3.0 Product Files may NOT work** for you.

For custom arrays:
- You must use **GenomeStudio** to build a **custom cluster file**
- You will receive:
  - A **custom manifest file**
  - Your **IDAT files**
  - A **sample sheet**

These are typically provided to you by the **array service provider** along with your data.

### Where standard GSA v3.0 product files live

The standard (non-custom) Infinium Global Screening Array v3.0 product files can be found here:

https://support.illumina.com/array/downloads.html

⚠️ These files are **only valid for standard arrays**, not custom designs.

---

## 5. GenomeStudio for custom cluster generation

To generate custom cluster files, you need **GenomeStudio**.

### GenomeStudio software
- Use **GenomeStudio Software 2.0.5**
- Download installer and release notes from Illumina Support

### Tutorial: creating custom cluster files

Illumina provides a step-by-step video tutorial:

**GenomeStudio™ Genotyping: Creating Custom Cluster Files for Infinium™ Arrays**  
[Watch the video](https://youtu.be/4JTrbMUbVN0?si=La2wOuZ5ypzjlsf-)

This process is required **before** downstream steps such as:
- gtc generation
- gtc → VCF conversion
- DRAGEN processing

---

## 6. Recommended reference for GenomeBuild 37

For Illumina genotyping arrays using **GenomeBuild 37**, use the Illumina-distributed GRCh37 reference:

```sh
ref_url=https://webdata.illumina.com/downloads/productfiles/microarray-analytics-array/GRCh37_genome.zip
wget $ref_url
unzip GRCh37_genome.zip
```

## 7. Why genome build consistency matters

Illumina tools (gtc → VCF, DRAGEN) are **strict** about:
- Chromosome naming
- Contig presence
- Coordinate system

A mismatched reference can cause:
- REF/ALT allele mismatches
- Missing contig errors
- Tool crashes or silent failures

Genome build consistency is **mandatory**, not optional.

---

## 8. Reference validation checklist (do not skip)

### 8.1 Integrity check (corruption)

If the FASTA is gzipped:

```bash
gunzip -t your_reference.fa.gz

If uncompressed, continue.

### 8.2 Index the reference (Required by downstream tools)

```sh
samtools faidx your_reference.fa
bwa index your_reference.fa
```

### 8.3 Inspect contig naming (critical)

```sh
awk '/^>/{print $1}' your_reference.fa | sed 's/^>//' | head -50
```

Check for:
- chr1 vs 1
- chrM vs MT
- Presence of ALT or random contigs

The contig set must match what your genotyping tools expect.

### 8.4 Quick sanity checks (Unexpected contig counts often indicate the wrong reference)

```sh
# number of contigs
cut -f1 your_reference.fa.fai | wc -l

# contig sizes
head your_reference.fa.fai

```
