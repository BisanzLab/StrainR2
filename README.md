# Background

Traditional methods for quantifying strain abundances in a microbiome, such as 16S rRNA sequencing, lack the resolution to differentiate strains and are limited to generalizing species. Metagenome assembled genomes have seen a recent explosion in popularity and can quantify abundance using a metric called FPKM: mapped fragments per thousand bases per million reads. The downfall of using FPKM is a bias from highly similar genomes getting mapped to less. 

Metagenomic methods that have been developed to remove the bias are inaccurate and scale poorly to larger communities. Abundances are off by magnitudes when strains have low coverage. They also produce these results after hours or days of computation with demanding resource usage. 

StrainR2 is a solution that quantifies strain abundances to within a couple of percent of the true value in a fraction of the time that traditional methods use. 


# Installation

StrainR2 is made for linux based operating systems. Dependencies are listed at the bottom of this README, but using conda is highly recommended for reproducibility.

### Option A: Bioconda (recommended)
Conda can be installed from the command line and is highly recommended for managing package dependencies.

```
conda create -n strainr2 -c bioconda -c conda-forge strainr2
conda activate strainr2
```

### Option B: Installation from source
To install the source code into a directory onto your computer, clone the source git repository:
```
git clone https://github.com/BisanzLab/StrainR2.git
```

Dependencies need to be installed according to versions listed at the bottom of this document. A .yml file provided in the git repository can be used to create an environment from scratch. 
```
conda env create -f StrainR2/strainr2.yml
conda activate strainr2
```
In addition subcontig.c will need to be compiled and executeable files will need executeable permissions:
```
gcc StrainR2/subcontig.c -o StrainR2/subcontig
chmod +x StrainR2/subcontig
chmod +x StrainR2/PreProcessR
chmod +x StrainR2/hashcounter.py
chmod +x StrainR2/StrainR
chmod +x StrainR2/Plot.R
```
Note that if StrainR2 is installed from source, the directory containing source files (the cloned git directory) needs to be added to your path every time a new shell is started. Replace `~/absolute/path/to` with the absolute path to the StrainR2 clone directory.
```
export PATH="$PATH:~/absolute/path/to/StrainR2"
```
<p>&nbsp;</p>

# Usage

StrainR2 is split into two steps: `PreProcessR` and `StrainR`. 

### PreProcessR 
`PreProcessR` creates a database for future runs of `StrainR` to use. It will split genome contigs into subcontigs to ensure similar genome build qualities, with contigs below a preset value (10kbp by default) being marked as "exluded". Excluded subcontigs have their k-mers marked as unique, but will not have their FUKM calculated. K-mers in eaach subconig are then compared to all other k-mers in the set of genomes to determine which ones are unique. Information about the count of unique k-mers is stored to KmerContent.report. In addition preprocessing will generate a `BBMap` index for later use. All of these files will be in the specified output directory.

`PreProcessR` will require the most memory. However, it only needs to run once for a given community of genomes and its output can be reused any number of times by `StrainR`. For large run sizes, it may be necessary to increase the size of swap space to facilitate memory needs. This would only be needed if `PreProcessR` crashes due to exceeding memory constraints.

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

Minimum size for subcontigs. Only contigs under the exclude size will be excluded. Default is 10000

**-s or --subcontigsize:**

This option overrides the default usage of the N50 build quality statistic, which StrainR2 calculates itself. StrainR2 will optimally create subcontigs with the given size as a maximum. Leaving the max size at N50 is recommended unless if the goal is to compare how changing max subcontig size affects the output. Smaller subcontigs make StrainR2 less sensitive to low abundance strains whereas larger subcontigs make StrainR2 too sensitive and it may report false abundances for absent strains due to index hopping and other sources of false reads. 

**-r or --readsize:**

The size of one end of your paired end reads. 150 by default.

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

**-o or --outdir:**

path to your output directory to contain normalized abundances. Default = current directory

**-p or --prefix:**

Name of community (used for output file names). Default = "sample"

**-t or --threads number:**

number of threads to use when running `fastp`, `BBMap`, and `samtools`. Maximum is 16, default = 8.

**-m or --mem number:**

gigabytes of memory to use when running `BBMap`. Default = 8.

**-s or --scratch string:** 

name of temporary directory used during run time [Default = tmp]

<p>&nbsp;</p>

# Outputs

### KmerContent.Report (output from PreProcessR) is formatted into the following columns:

SubcontigID: Formatted as follows: \<FILE_NAME\>;\<CONTIG_HEADER_NAME\>;\<START-NUCLEOTIDE_STOP-NUCLEOTIDE\>;\<SUBCONTIG_LENGTH\>

StrainID: File name for strain

ContigID: Contig header from .fasta or .fna file

Start_Stop: Nucleotides at which the subcontig starts and stops. Counting starts from the start of a strain genome file to the end.

Length: Number of basepairs in subcontig

Nunique: Number of unique kmers in subcontig


### The .abundance file is formatted into the following columns:

StrainID, ContigID, Start_Stop, Unique_Kmers, Length_Contig: Same as KmerContent.report file

Bases, Coverage, Mapped_Reads, Mapped_Frags, Total_Mapped_Reads_In_Sample: Outputs from `BBMap` in .rpkm file

FUKM: StrainR2-calculated strain abundance

### The abundance_summary.tsv file is formatted into the following columns:

StrainID: Name of the fasta file for which summary statistics will be provided.

mean_FUKM, median_FUKM, sd_FUKM: Mean, median, and standard deviation of all subcontig FUKMs in the strain

subcontigs_detected: The number of subcontigs that had an FUKM>0

subcontigs_total: The total number of subcontigs for the strain.

<p>&nbsp;</p>

In addition, `StrainR` provides a plot for FUKM abundances. Median FUKM is the recommended measure of strain abundance.


<p>&nbsp;</p>

# Example Run

The following example of running StrainR2 on mock reads of the sFMT community serves as an example of usage as well as a way to verify reproducibility. Genomes of sFMT strains as well as mock read files are provided in `example_data`. The following command creates a StrainR2 data base which can be used by the `StrainR` command for all subsequent read analysis. 
```
PreProcessR -i StrainR2/example_data/sFMT/ -o <output/directory/of/your/choosing>
```
This will generate a StrainR2 database with default options (N50 subcontig size, 10kbp exclude size, 150bp reads), and put it in the output directory specified by the `-o` option.

`StrainR` can be run multiple times on the output of `PreProcessR`. An example using mock reads where each strain has equal coverage is given below. The data used for runs as well as the corresponding output can be found in `example_data`.

```
StrainR -1 StrainR2/example_data/mock_reads/shallow_uniform_R1.fastq.gz -2 StrainR2/example_data/mock_reads/shallow_uniform_R2.fastq.gz -r <your/preprocessr/output/directory> -p shallow_uniform -o <output/directory/of/your/choosing>
```
Abundance data will be put in the directory specified after `-o`. This includes a plot of abundances, which in the case of the provided mock reads will be very close to uniform. Due to the provided reads having low coverage to keep file sizes small, subcontig FUKMs will have a larger spread than usual. We recommend using median_FUKM as the measure of abundance for strains.

<p>&nbsp;</p>

# Dependencies

 * sourmash >4.0.0
 * BBMap
 * fastp
 * samtools
 * R
 * R tidyverse

# Credits:
Paper: \<Link to paper to be added when published\>

For any issues or questions contact kerim@heber.org

<p>&nbsp;</p>

