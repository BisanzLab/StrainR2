[![Testing](https://github.com/BisanzLab/StrainR2/actions/workflows/testing.yml/badge.svg)](https://github.com/BisanzLab/StrainR2/actions/workflows/testing.yml)
[![StrainR2 Version](https://anaconda.org/bioconda/strainr2/badges/version.svg)](https://anaconda.org/bioconda/strainr2)
[![Downloads](https://anaconda.org/bioconda/strainr2/badges/downloads.svg)](https://anaconda.org/bioconda/strainr2)
# Background

Traditional methods for quantifying strain abundances in a microbiome, such as 16S rRNA sequencing, lack the resolution to differentiate strains and are limited to generalizing species. Shotgun metagenomic sequencing offers an alternative, but unnormalized abundances such as FPKM have a bias from similar genomes getting fewer unique mappings. 

Other metagenomic methods developed to remove this bias, such as NinjaMap or StrainR1, are inaccurate and scale poorly to larger communities. Abundances are off by magnitudes especially when strains have low coverage, and computational resources scale poorly.

StrainR2 is a solution that quantifies strain abundances using shotgun metagenomic sequencing to within a couple of percent of the true value in a fraction of the time that traditional methods use. 


# Installation

StrainR2 is made for unix-like operating systems. Dependencies are listed at the bottom of this README, but using conda is highly recommended for reproducibility.

### Option A: Bioconda (recommended)
Conda can be installed from the command line and is highly recommended for managing package dependencies.

```
conda create -n strainr2 -c bioconda -c conda-forge strainr2
conda activate strainr2
```

### Option B: Docker
Docker is another package manager that can be used to ensure StrainR2 can be run in a replicable manner. 
```
docker pull quay.io/biocontainers/strainr2:<tag>
```
Where \<tag\> needs to be replaced by any valid tag ([see valid StrainR2 tags](https://quay.io/repository/biocontainers/strainr2?tab=tags)).

### Option C: Installation from source
To install the source code into a directory onto your computer, clone the source git repository:
```
git clone https://github.com/BisanzLab/StrainR2.git
```

Dependencies need to be installed according to versions listed at the bottom of this document.

Files can be compiled using make
```
cd StrainR2/src
make release
```

Note that if StrainR2 is installed from source, code needs to be run by referencing the appropriate files in the `src` directory directly.

<p>&nbsp;</p>

# Usage

StrainR2 is split into two steps: `PreProcessR` and `StrainR`. 

### PreProcessR 
`PreProcessR` creates a database for future runs of `StrainR` to use. It will split genome contigs into subcontigs to ensure similar genome build qualities, with the contigs that are below a preset value (10kbp by default) being marked as "exluded". Information about the count of unique k-mers is stored to KmerContent.report. In addition preprocessing will generate a `BBMap` index for later use. All of these files will be in the specified output directory.

`PreProcessR` only needs to run once for a given community of genomes and its output can be reused any number of times by `StrainR`. For large run sizes, it may be necessary to increase the size of swap space to facilitate memory needs. This would only be needed if `PreProcessR` crashes due to exceeding memory constraints.

The `PreProcessR` command can be invoked from the command line as follows:
```
PreProcessR -i <PATH_TO_GENOME_FILES> [OPTIONS]
```

#### Options explanation:

**-i or --indir [REQUIRED]:**

The path to a directory containing only the genome files to be quantified. Accepted file formats are .fna and .fasta

**-o or --outidr:**

Path to desired directory for output files to be stored in. Directory does not need to be made before running `PreProcessR`. If the output directory already exists before a run, make sure it is empty. Default is `./StrainR2DB`

**-e or --excludesize:**

Minimum size for subcontigs. Only contigs under the exclude size will be excluded. Default = 10,000

**-s or --subcontigsize:**

This option overrides the default usage of the minimum N50 build quality statistic of all inputted genomes, which StrainR2 calculates itself. StrainR2 will optimally create subcontigs with the given size as a maximum. Leaving the max size as default is recommended. Smaller subcontigs make StrainR2 less sensitive to low abundance strains whereas larger subcontigs make StrainR2 too sensitive and it may report false abundances for absent strains due to index hopping, read errors, and other sources of false reads. 

**-r or --readsize:**

The size of one end of your paired end reads. 150 by default. To be used to calculate a k-mer size.

<p>&nbsp;</p>


### StrainR
`StrainR` takes paired end reads and normalizes the abundance of strains using the output prepared in the `PreProcessR` step. It will generate a .sam, .rpkm, and .abundance file. The .sam and .rpkm files are generated using `BBMap` and the .abundance file normalizes .rpkm output by the data in the previously generated KmerContent.report file. In addition, `StrainR` will generate a plot of abundance.

The `StrainR` command can be invoked from the command line as follows:
```
StrainR -1 <PATH_TO_FORWARD_READS> -2 <PATH_TO_REVERSE_READS> -r <PATH_TO_OUTPUT_OF_PREPROCESSR>[OPTIONS]
```

#### Options explanation:
**-1 or --forward [REQUIRED]:**

path to forward reads (reads do not have to be in .gz)

**-2 or --reverse [REQUIRED]:**

path to reverse reads (reads do not have to be in .gz)

**-r or --reference [REQUIRED]:**

path to the output directory generated by `PreProcessR` (the directory that was the outout for `PreProcessR`)

**-c or --weightedpercentile:**

The weighted percentile of a strain's subcontig's FUKMs to be used for final abundance calculation (weighted by the number of unique k-mers). Adjusting this to be higher will increase sensitivity but make false positives more likely. Conversely, decreasing this number will increase the likelihood of false negatives. Default = 60 

**-s or --subcontigfilter:**

The percentage of subcontigs that should get filtered out on the basis of having too few unique k-mers. Default = 0

**-o or --outdir:**

path to your output directory to contain normalized abundances. Default = current directory

**-p or --prefix:**

Name of community (used for output file names). Default = "sample"

**-t or --threads number:**

number of threads to use when running `fastp`, `BBMap`, and `samtools`. Maximum is 16, default = 8

**-m or --mem number:**

gigabytes of memory to use when running `BBMap`. Default = 8

<p>&nbsp;</p>

# Outputs

### KmerContent.Report (output from PreProcessR) is formatted into the following columns:

SubcontigID: Formatted as follows: \<FILE_NAME\>;\<CONTIG_HEADER_NAME\>;\<START-NUCLEOTIDE_STOP-NUCLEOTIDE\>;\<SUBCONTIG_LENGTH\>

StrainID: File name for strain

ContigID: Contig header from .fasta or .fna file

Start_Stop: Nucleotides at which the subcontig starts and stops. Counting starts from the start of a strain genome file to the end.

Length: Number of basepairs in subcontig

Nunique: Number of unique k-mers in subcontig


### The .abundance file is formatted into the following columns:

StrainID, ContigID, Start_Stop, Unique_Kmers, Length_Contig: Same as KmerContent.report file

Bases, Coverage, Mapped_Reads, Mapped_Frags, Total_Mapped_Reads_In_Sample: Outputs from `BBMap` in .rpkm file

FUKM: StrainR2-calculated strain abundance

### The abundance_summary.tsv file is formatted into the following columns:

StrainID: Name of the fasta file for which summary statistics will be provided.

weighted_percentile_FUKM, median_FUKM, sd_FUKM: Weighted percentile (per `StrainR`'s `-c` option), median, and standard deviation of all subcontig FUKMs in the strain

subcontigs_detected: The number of subcontigs that had an FUKM>0

subcontigs_total: The total number of subcontigs for the strain

percent_detected: subcontigs_detected / subcontigs_total

percent_abundance: weighted_percentile_FUKM / sum of all weighted_percentile_FUKM in the community

<p>&nbsp;</p>

In addition, `StrainR` provides a plot for FUKM abundances. Weighted percentile FUKM is the recommended measure of strain abundance.

<p>&nbsp;</p>

# Dependencies

 * BBMap
 * fastp
 * GNU make (if compiling from source)
 * samtools
 * R
 * R optparse
 * R tidyverse
 * zlib

# Credits:
Paper: \<Link to paper to be added when published\>

For any issues or questions contact kerim@heber.org

<p>&nbsp;</p>

