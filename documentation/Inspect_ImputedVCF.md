

Column by column (what matters):

10 → chromosome

6088743 → position (GRCh37)

rs12722508 → rsID (this is what your script reads as rsID)

A / T → reference / alternate allele

PASS → passed filters

INFO:

AF=0.0995 → alt allele frequency (~10%)

DR2=1 → high imputation quality

IMP → imputed variant

GENE=IL2RA → annotated gene

FORMAT = GT:DS

GT → genotype (phased)

0|0 = hom ref

0|1 or 1|0 = heterozygous

1|1 = hom alt

DS → dosage (expected ALT allele count)

0, 1, 2 (or decimals for imputed calls)