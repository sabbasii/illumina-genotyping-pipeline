#

# Initialize the R session
rm(list = ls())

suppressPackageStartupMessages({
	library(tidyverse)
	library(conflicted)
})

conflict_prefer_all("dplyr", quiet = TRUE)

args <- commandArgs(trailingOnly = TRUE)
report_file <- args[1]

# Load in target report
targets <- read_tsv(
	file = report_file,
	col_types = cols(
		REF = col_character(),
		ALT = col_character()
	)
)

#
cat(
	"\n",
	"Cleaning report:...",
	"\n"
)

# Clean genotype data by removing unnecessary columns and renaming for clarity
targets_cleaned <- targets |>
	select(
		!c(QUAL, FILTER, FORMAT)
	) |>
	rename(
		CHR = "#CHROM"
	)

# Create vector of sample names
samples <- targets_cleaned |>
	select(
		starts_with("UASG")
	) |>
	colnames()

# Record total number of markers for reference
total_markers <- nrow(targets_cleaned)

# Split the sample values and return only the genotypes
targets_parsed <- targets_cleaned |>
	mutate(
		across(
			all_of(samples),
			~ str_split_i(.x, ":", 1)
		)
	)

# Separate sample markers from imputed.
measured <- targets_parsed |>
	filter(
		! str_detect(INFO, ";IMP")
	) |>
	separate_wider_delim(
		cols = INFO,
		delim = ";",
		names = c("AC", "AN", "GENE"),
		too_few = ("align_start")
	) |>
	mutate(
		AC = as.numeric(str_replace(AC, "AC=", "")),
		AN = as.numeric(str_replace(AN, "AN=", "")),
		GENE = str_replace(GENE, "GENE=", ""),
		AF = round((AC / AN), 2),
		DR2 = NA,
		.after = ALT,
	) |>
	select(
		! c(AC, AN)
	)

total_measured <- nrow(measured)

imputed <- targets_parsed |>
	filter(
		str_detect(INFO, ";IMP")
	) |>
	separate_wider_delim(
		cols = INFO,
		delim = ";",
		names = c("DR2", "AF", "IMP", "GENE"),
		too_few = ("align_start")
	) |>
	mutate(
		DR2 = as.numeric(str_replace(DR2, "DR2=", "")),
		AF = round(as.numeric(str_replace(AF, "AF=", "")), 2),
		GENE = str_replace(GENE, "GENE=", "")
	) |>
	select(
		! IMP
	)

total_imputed <- nrow(imputed)

#
cat(
	"\n",
	"Summary of markers in report:",
	"\n",
	"Total Markers: ",
	total_markers,
	"; Measured Markers: ", 
	total_measured, 
	"; Imputed Markers: ", 
	total_imputed,
	"\n"
)

# Combine and sort the dataframes, for multi-allelics, put the higher AF first.
targets_combined <- measured |>
	bind_rows(
		imputed
	) |>
	arrange(
		CHR,
		POS,
		desc(AF)
	)

# Convert the VCF notation of 0|0, 1|0, etc. into number of ALT alleles.
allele_counts <- targets_combined |>
	mutate(
		across(
			all_of(samples),
			~ str_count(.x, "1")
		)
	)

write_tsv(
	allele_counts,
	file = report_file
)
