## 1. Introduction

Genotype imputation is a standard step after initial genotyping QC. The goal is to infer _missing_ or _untyped variants_ using a reference population, while keeping the original measured genotypes unchanged.

This document describes **what imputation is**, **why we do it**, and **how it is handled in this pipeline**. It focuses on practical steps, file requirements, and how to inspect the output.

The workflow here assumes:
- Genotypes have already passed basic QC and normalization
- Data are aligned to a known genome build (e.g. GRCh37)
- Only autosomal and X chromosome variants are considered for imputation

---

## 2. Scope and goals of imputation

**Scope**
- Applied after genotype QC and normalization
- Uses external reference panels to infer missing or untyped variants
- Limited to autosomal and X chromosomes
- Does not modify original measured genotypes

**Goals**
- Increase variant density for downstream analysis
- Fill in missing genotypes where direct calls are unavailable
- Maintain compatibility with reference-based analyses (GWAS, PRS, eQTL)
- Preserve a clear distinction between measured and imputed variants

**Out of scope**
- Structural variants
- Y and mitochondrial chromosomes
- Pseudoautosomal regions (PAR)
- Variant discovery or re-calling of genotypes

---

## 3. What imputation is

Imputation estimates genotypes at variants that were not directly measured on the genotyping array.

This is done by comparing observed haplotypes in the study samples to haplotypes in a reference population. If a pattern in the study data matches a pattern in the reference, missing variants along that pattern can be inferred.

The output genotypes are probabilistic, not direct measurements. Each imputed variant is accompanied by:
- An estimated allele dosage
- A quality metric describing confidence in the estimate

Imputation does not create new biological signal. It leverages linkage disequilibrium between nearby variants to infer likely genotypes.


**Reference**
- Browning BL, Zhou Y, Browning SR. *A one-penny imputed genome from next-generation reference panels.* Am J Hum Genet. 2018.  
- https://faculty.washington.edu/browning/beagle/beagle.html

---

## 4. Data and reference prerequisites

Imputation relies on matching **haplotypes** between the study data and a reference panel.  

<div style="background-color: #e5cad0; color:#000000; padding:12px 14px; border-left:4px solid #8e5b87; border-radius:6px; margin:12px 0;">

A *haplotype* is the ordered set of alleles carried together on a single chromosome copy.  
Because genotyping arrays measure variants without knowing which chromosome copy they came from, haplotypes must first be inferred through **phasing**.
</div>

<div style="background-color: #e5cad0; color:#000000; padding:12px 14px; border-left:4px solid #8e5b87; border-radius:6px; margin:12px 0;">

**_Phasing_** is the process of assigning alleles to the maternal or paternal chromosome, producing two haplotypes per individual. Imputation operates on these phased haplotypes, not on unphased genotypes.
</div>

---

### Target genotype data

The target dataset is the study genotype file produced by the array pipeline.

It must:
- Be generated **after genotype QC and normalization**
- Use a **single genome build** (e.g. GRCh37)
- Be stored as a **compressed and indexed VCF or BCF**

In this workflow, the target input is:  
`genotypes_normalized.bcf`


This file provides the observed variants used to:
- Phase the samples
- Match sample haplotypes to reference haplotypes

Genotypes in this file are treated as fixed observations and are not re-called during imputation.

Before imputation, the target data must be compatible with the reference in terms of:
- Chromosome names (exact string match)
- Genome build
- Supported chromosomes (autosomes and X only)
- Correct ploidy representation for chromosome X

Failure to meet any of these conditions will cause BEAGLE to error or silently drop data.

---

### Reference panel

The reference panel provides **phased haplotypes** derived from a large external population.

The reference must:
- Be phased
- Match the **same genome build** as the target data
- Use the **same chromosome naming scheme**
- Be compressed and indexed

In this workflow, the reference comes from the `BEAGLE 1000 Genomes Phase 3 dataset`, distributed as per-chromosome VCF files.   
These are concatenated into a single reference file:  
`concatenated_ref.vcf.gz`


Only variants present in the reference panel can be imputed. Variants present in the target dataset but absent from the reference are temporarily removed during imputation and handled later.

---

### Genetic map files

Genetic map files describe recombination rates along each chromosome. These rates guide phasing and imputation by defining how likely it is for haplotypes to break and recombine between variants.

The map files must:
- Match the genome build (e.g. GRCh37)
- Cover all chromosomes included in the analysis
- Be provided in PLINK `.map` format

In this workflow, per-chromosome map files are downloaded from the BEAGLE site and concatenated into:  
`concatenated_map_file.map`


Pseudoautosomal regions are excluded and not imputed.

---

### Variant overlap requirement

BEAGLE requires **variant overlap** between the target dataset and the reference panel to align haplotypes.

- Shared variants are used as anchors to match haplotypes
- Reference-only variants are **imputed into the target samples**
- Target-only variants are **not carried through BEAGLE output**

This is why the imputed dataset contains:
- New imputed variants from the reference
- A reduced set of original variants

Measured variants dropped at this stage are later merged back during dataset reconstitution.  

---

**Reference**
- Browning BL, Browning SR. *Genotype Imputation with Millions of Reference Samples.*  
  American Journal of Human Genetics, 2016.  
  https://faculty.washington.edu/browning/beagle/beagle.html