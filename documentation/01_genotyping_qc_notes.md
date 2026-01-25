# Illumina Cluster QC Notes

**Author:** Sima Abbasi-Habashi

---

## Table of Contents

1. Introduction  
2. Genotype Clusters (AA, AB, BB)  
3. Cluster Quality Scores  
   3.1 GenTrain_Score  
   3.2 Orig_Score  
   3.3 Cluster_Sep  
   3.4 Edited Flag  
4. Training Set Genotype Counts (N_AA, N_AB, N_BB)  
5. Cluster Signal Metrics  
   5.1 R (Total Intensity)  
   5.2 THETA (Allele Balance)  
6. Cluster Statistics per Genotype  
   6.1 meanR and devR  
   6.2 meanTHETA and devTHETA  
7. Genotype Calling Geometry (R vs THETA)  
8. Log R Ratio (LRR)  
9. Allele Frequency (AF)  
10. Call Rate and Missingness  
11. Hardy–Weinberg Equilibrium (HWE)  
12. Sample-Level QC  
   12.1 Sample Call Rate  
   12.2 LRR Standard Deviation  
   12.3 Heterozygosity Rate  
13. SNP-Level QC Thresholds  
   13.1 Call Rate  
   13.2 Minor Allele Frequency (MAF)  
   13.3 Hardy–Weinberg Equilibrium  
   13.4 Cluster-Based Filters  

---

## 1. Introduction

This document summarizes the key concepts and metrics used to assess _SNP clustering and genotype quality_ in Illumina genotyping arrays.  
It focuses on

+ how allele signal intensities are translated into genotype calls
+ how cluster quality is evaluated
+ which statistics are used to identify reliable versus problematic SNPs.

---

## 2. Genotype Clusters (AA / AB / BB)

For each SNP, Illumina arrays measure fluorescence intensities for two alleles:

- A_intensity
- B_intensity

Each sample is represented as a point in intensity space based on these two values. When all samples are plotted for a SNP, they form three groups corresponding to the three possible genotypes:

- __AA__: high A_intensity, low B_intensity
- __AB__: intermediate A_intensity and B_intensity
- __BB__: low A_intensity, high B_intensity

Each group of points is a **genotype cluster**. A high-quality SNP shows three tight, well-separated clusters. Poor-quality SNPs may show overlapping clusters, elongated clusters, missing clusters, or excessive scatter, which leads to unreliable genotype calls.

---

## 3. Cluster Quality Scores

Illumina provides several metrics to quantify the quality of genotype clustering for each SNP. These scores summarize cluster separation, tightness, and overall reliability.

### 3.1 GenTrain_Score

GenTrain_Score is a value between 0 and 1 that reflects overall cluster quality as assessed by Illumina’s GenTrain algorithm. It captures how tight the clusters are and how well they are separated from each other. Higher values indicate more reliable genotype calls.

Typical interpretation:

- 0.7–0.8 and above: good-quality SNP
- Below ~0.3: unreliable clustering

### 3.2 Orig_Score

Orig_Score is the GenTrain score before any manual cluster editing. It allows comparison between the original automated clustering and the final clustering. A higher final GenTrain_Score than Orig_Score indicates that manual intervention improved cluster quality.

### 3.3 Cluster_Sep

Cluster_Sep measures how far apart the genotype clusters are, on a scale from 0 to 1. This metric focuses on separation rather than cluster tightness. Higher values indicate minimal overlap between AA, AB, and BB clusters, while lower values indicate blurred or overlapping clusters.

### 3.4 Edited Flag

The Edited flag indicates whether the SNP’s clusters were manually adjusted during curation. This does not automatically imply poor quality, but edited SNPs should be reviewed with more caution, especially if other quality metrics are borderline.

---

## 4. Training Set Genotype Counts (N_AA, N_AB, N_BB)

`N_AA`, `N_AB`, and `N_BB` represent the number of samples used to define the `AA`, `AB`, and `BB` genotype clusters during cluster training.  
These counts reflect _how well each genotype class is represented_ when the cluster model is built.

`Low counts` in any genotype reduce cluster stability and increase uncertainty in genotype boundaries.

+ Extremely unbalanced counts, such as very few heterozygotes or a missing genotype class, may indicate a rare variant, a poorly performing SNP, or systematic calling issues.

`Training set genotype counts` should be interpreted alongside cluster quality scores and visual inspection of cluster plots to assess overall reliability.

---

## 5. Cluster Signal Metrics

Cluster signal metrics describe the underlying fluorescence intensity structure used for genotype calling. Illumina primarily relies on two signal dimensions: total intensity and allele balance.

### 5.1 R (Total Intensity)

`R` represents the total signal intensity measured at a SNP and is defined as the sum of the two allele intensities:

`R = A_intensity + B_intensity`

`R` reflects the overall strength of the hybridization signal.  
`Low R values` indicate weak signal and are often associated with noisy or unreliable genotype calls.  
`High R values` indicate strong, reliable measurements.  
Quality control procedures commonly filter SNPs or samples with consistently low R.

### 5.2 THETA (Allele Balance)

`THETA` represents the relative balance between the two allele intensities and _ranges from 0 to 1_.  
It captures how signal is distributed between alleles:

- AA genotypes have THETA values near 0
- AB genotypes have THETA values near 0.5
- BB genotypes have THETA values near 1

`THETA` is the primary dimension used to distinguish genotypes, while `R` mainly reflects confidence in the measurement.

---

## 6. Cluster Statistics per Genotype

For each genotype cluster (AA, AB, BB), Illumina reports summary statistics that describe the location and spread of the cluster in signal space.   These metrics quantify cluster stability and consistency across samples.

### 6.1 meanR and devR

`meanR` is the average total intensity (R) of samples within a genotype cluster.  
`devR` describes the spread of R values around this mean.

`Low devR` indicates a _tight cluster_ with consistent signal strength, while ``high devR` indicates variable intensity and reduced confidence in genotype calls.

### 6.2 meanTHETA and devTHETA

`meanTHETA` is the average allele balance for a genotype cluster.  
`devTHETA` measures the spread of THETA values within the cluster.

`Low devTHETA` indicates a well-defined genotype cluster, while `high devTHETA` indicates _fuzzy or overlapping clusters_.

* Tight THETA distributions are critical for accurate genotype separation.

---

## 7. Genotype Calling Geometry (R vs THETA)

Genotype calling on Illumina arrays is based on the joint geometry of R and THETA. THETA determines the genotype assignment, while R determines the confidence of that assignment.

Clusters are primarily separated along the THETA axis, which distinguishes AA, AB, and BB genotypes based on allele balance. The R axis reflects signal strength and does not define genotype boundaries, but influences how reliable a call is.

High-confidence genotype calls have THETA values close to the expected cluster centers and sufficiently high R values. Calls with appropriate THETA but low R are more prone to noise and misclassification. As a result, quality control focuses on THETA for genotype discrimination and R for call reliability.

---

## 8. Log R Ratio (LRR)

`Log R Ratio (LRR)` is a normalized measure of total signal intensity used primarily for copy-number analysis rather than genotype calling.  
It compares the `observed intensity` at a SNP to the `expected intensity` derived from Illumina’s reference model.

LRR is defined as:

`LRR = log2(R_observed) − log2(R_expected)`

`R_observed` is the measured total intensity (A_intensity + B_intensity).  
`R_expected` is obtained from reference cluster models stored in the manifest or cluster (EGT) files and is not estimated from the study data.

Interpretation:

- LRR ≈ 0 indicates normal diploid copy number
- LRR < 0 indicates reduced intensity, consistent with deletions
- LRR > 0 indicates increased intensity, consistent with duplications

`LRR` is an intensity-level metric and is not included in standard genotyping VCF or BCF files. It must be exported from _GenomeStudio_ or generated through copy-number analysis pipelines.

---

## 9. Allele Frequency (AF)

`Allele frequency (AF)` describes how common the alternate allele is in a population.  
It is calculated as the proportion of all observed alleles that are the alternate allele.

AF is defined as:

`AF = number of alternate alleles / total number of alleles`

Interpretation:

- AF near 0 indicates a `rare variant`
- AF near 0.5 indicates a `common variant`
- AF near 1 indicates the alternate allele is the `major allele`

`AF` is a population-level metric and is used for variant characterization and quality control, rather than for assessing individual genotype confidence.

---

## 10. Call Rate and Missingness

`Call rate` and `missingness` quantify how often genotypes are successfully assigned, at both the `SNP level` and the `sample level`.  
They are core quality control metrics in array-based genotyping.

### Definitions

- **Call rate**: proportion of non-missing genotype calls
- **Missingness**: proportion of missing genotype calls

These metrics are complementary:

`Missingness = 1 − Call rate`

### How to Calculate

**Per SNP:**

`Call rate (SNP)` =  number of samples with a genotype call / total number of samples

`Missingness (SNP)` =  number of missing genotypes / total number of samples

**Per sample:**

`Call rate (sample)` =  number of SNPs successfully genotyped / total number of SNPs

`Missingness (sample)` =  number of missing SNPs / total number of SNPs

### Where to Get It (Which Files)

- **VCF / BCF (post-genotyping)**  
   Missing genotypes are encoded as `./.`  
   Call rate and missingness can be computed directly from these files and are commonly summarized by downstream tools.
- __PLINK files (PGEN / BED)__  
   Call rate and missingness are typically computed during _QC using PLINK_, which reports:

   - SNP-level missingness
   - Sample-level missingness

- __GenomeStudio reports (pre-VCF)__  
   Call rate can also be exported directly from _GenomeStudio_ at the SNP or sample level, prior to VCF/BCF generation.

### Note:

- `Low SNP call rate` suggests poor clustering, weak signal, or technical failure.
- `Low sample call rate` suggests poor DNA quality or assay failure.
- Call rate filters are usually applied before downstream analyses such as association testing or imputation.

---

## 11. Hardy–Weinberg Equilibrium (HWE)

Hardy–Weinberg equilibrium (HWE) assesses whether genotype frequencies at a SNP are consistent with expectations under random mating. It is used as a quality control check to identify genotyping errors or problematic variants.

### Concept

For a bi-allelic SNP with allele frequencies p (REF) and q (ALT), where p + q = 1, expected genotype frequencies are:

- AA: p²
- AB: 2pq
- BB: q²

Observed genotype counts are compared to these expectations.

### How It Is Evaluated

HWE is typically tested using a statistical test (chi-square or exact test) that produces a p-value. A low p-value indicates deviation from equilibrium.

### Where to Get It (Which Files)

- **VCF / BCF**  
   HWE is not stored directly but is computed from genotype calls.
- **PLINK files (PGEN / BED)**  
   HWE statistics are routinely calculated during QC and reported per SNP.

### Practical Use in QC

- HWE is usually evaluated in controls or unaffected samples.
- Strong deviations may indicate genotyping errors, poor clustering, or allele-specific bias.
- Genuine biological effects can also cause deviation, so HWE failures should be interpreted in context.

### Common Practice

- Apply HWE filtering as part of SNP-level QC.
- Use stricter thresholds in controls and more lenient thresholds in cases, depending on study design.

---

## 12. Sample-Level QC (Call Rate, LRR SD, Heterozygosity)

`Sample-level QC` evaluates _overall data quality per individual_ and is used to identify:

+ samples with poor DNA quality
+ technical failures
+ abnormal genomic patterns

### 12.1 Sample Call Rate

Sample call rate measures the proportion of SNPs successfully genotyped for a given sample.

Formula:

`Sample call rate =  number of non-missing genotype calls / total number of SNPs`

`Missingness = 1 − call rate`

**Source files:**

- VCF / BCF (missing genotypes coded as `./.`)
- PLINK files (reported directly during QC)
- GenomeStudio sample reports (pre-VCF)

Low call rate samples are typically excluded early in QC.

### 12.2 LRR Standard Deviation (LRR SD)

`LRR SD` measures the variability of Log R Ratio across the genome for a sample.  
It reflects `signal noise`.

Formula (conceptual):

`LRR SD = standard deviation of LRR values across SNPs`

**Source files:**

- GenomeStudio LRR export
- CNV pipeline outputs

Interpretation:

- Low LRR SD → stable signal
- High LRR SD → noisy sample, potential hybridization or DNA quality issues

LRR SD is mainly used in copy-number and array QC, not genotype calling.

### 12.3 Heterozygosity Rate

`Heterozygosity rate` measures the proportion of heterozygous genotypes in a sample.

Formula:

`Heterozygosity rate =  number of heterozygous SNPs / number of non-missing SNPs`

**Source files:**

- VCF / BCF
- PLINK files (reported during QC)

Interpretation:

- Excessively high heterozygosity may indicate _contamination_
- Excessively low heterozygosity may indicate _inbreeding_ or _technical artifacts_

Samples with heterozygosity far from the cohort mean are typically flagged for review or removal.

---

## 13. SNP-Level QC Thresholds (Practical Cutoffs)

SNP-level QC applies filters _to remove variants with unreliable genotype calls or technical artifacts_.  
Thresholds may vary by study, but the following are commonly used in practice.

### 13.1. Call Rate / Missingness

- Call rate ≥ 0.95–0.98
- Missingness ≤ 0.02–0.05

Low call rate SNPs often reflect poor clustering or weak signal.

---

### 13.2. Minor Allele Frequency (MAF)

`Minor allele frequency (MAF)` is the frequency of the less common allele at a SNP in the study population.

`MAF = min(AF, 1 − AF)`

+ You always report the frequency of the `less common allele`, regardless of which allele is labeled REF or ALT.

`AF` = frequency of the alternate allele
`1 − AF` = frequency of the reference allele

+ The minor allele frequency (MAF) is the smaller of the two.
   Examples:
   AF = 0.10 → MAF = min(0.10, 0.90) = 0.10
   AF = 0.80 → MAF = min(0.80, 0.20) = 0.20
   AF = 0.50 → MAF = min(0.50, 0.50) = 0.50

Common cutoffs:

- MAF ≥ 0.01 for standard QC
- MAF ≥ 0.05 for stricter analyses

Very rare variants are more prone to unstable clustering and genotyping errors.

---

### 13.3. Hardy–Weinberg Equilibrium (HWE)

- Controls: p ≥ 1e−6 (common)
- Cases: p ≥ 1e−10 or no filter, depending on design

## Strong deviations may indicate genotyping problems rather than biology.

### 13.4. Cluster Quality Metrics (if available)

- GenTrain_Score ≥ 0.7
- Cluster_Sep high (no strict cutoff, but low values are flagged)
- Avoid SNPs with extreme imbalance in N_AA / N_AB / N_BB

---

### 13.5. Additional Common Filters

- Remove SNPs with ambiguous clustering (visual inspection)
- Remove SNPs with strand or annotation inconsistencies
- Exclude SNPs failing platform-specific QC flags

QC thresholds should be chosen consistently and documented clearly, as they directly affect downstream analyses.