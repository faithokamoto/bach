# bach
Biased Allele Contacts in Hi-C

BME 230A class project, heavily based on a [rotation project][ProjectSummary]
performed in [Dr. Richard "Ed" Green's lab][PGL] with Dr. Merly Escalona.

## Installation

### Get `bach`

Clone `bach` from GitHub:

```
git clone https://github.com/faithokamoto/bach
cd bach # Enter directory
chmod +x ./bach # Make executable
```

Test by:
```
./bach --help
```

### Dependencies

`bach` has the following dependencies, documented in `environment.yml`:
* Python 3 (tested on version `3.12.8`) [Installation instructions][PythonInstall]
* `pysam` (tested on version `0.23.0`) [Installation instructions][pysamInstall]
* `SciPy` (tested on version `1.15.2`) [Installation instructions][ScipyInstall]

> [!NOTE]
> Of these dependencies, the one most likely to have an issue is `SciPy`.
> Versions before `1.10.0` lack `scipy.stats.binomtest`. They instead use
> `scipy.stats.binom_test`, which has since been [deprecated][Deprecation].

Alternatively, if you have [conda][CondaInstall], simply run this from inside of
the `bach` directory:
```
conda env create --file environment.yml
conda activate bach
```

## Test data

To run `bach` on the provided test data, use 
```
./bach --max-opposite 0 --max-neutral 0 --max-drop 1 \
    --window-width 15 --move-step 5 \
    -v test/snps.vcf -d test/bams -o using_step.csv

./bach --max-opposite 0 --max-neutral 0 --max-drop 1 \
    --window-width 15 --round-to 5 \
    -v test/snps.vcf -d test/bams -o using_round.csv
```

The output files should correspond to the ones in `test`

To generate new test data, use 
```
python scripts/simulate_reads.py --depth 30 \
    --biased NC_000020.11:1182052 --biased NC_000020.11:1182049 \
    test/sequence.fasta test/snps.vcf test/bams
```

## How `bach` works

### Help message

```
usage: bach [-h] [-p #] [-v VCF] [-d BAM_DIR] [-o OUTPUT] [-e .<ext>] [-O 0.*]
            [-N 0.*] [-D #] [-w Mbp] [-s Mbp] [-r Mbp]

Find allele-biased segments. Compare mates of reads with SNP alleles.

options:
  -h, --help            show this help message and exit
  -p #, --print-progress #
                        Progress update frequency (# variants processed)

filepaths:
  -v VCF, --vcf VCF     Path to the SNP VCF
  -d BAM_DIR, --bam-dir BAM_DIR
                        Directory with BAM files
  -o OUTPUT, --output OUTPUT
                        Output CSV path
  -e .<ext>, --bam-ext .<ext>
                        Extension for BAMs (<sample>.ext)
                        [default=".pos.sorted.bam"]

bias parameters:
  -O 0.*, --max-opposite 0.*
                        Max share [0,0.5) of heterozygotes with opposite bias
                        [default=0]
  -N 0.*, --max-neutral 0.*
                        Max share [0,1) of heterozygotes with 0 change
                        [default=0]
  -D #, --max-drop #    Max count of homozygous/missing samples [default=0]

window parameters:
  -w Mbp, --window-width Mbp
                        Window width, in Mbp, to scan with (Consecutive windows
                        are merged)
  -s Mbp, --move-step Mbp
                        Mbp step size, in Mbp. Use one of -s or -r
  -r Mbp, --round-to Mbp
                        Mbp to round read positions to in empirical window
                        choice. Use one of -s or -r
```

### Workflow

#### Background

[Chromatin], the physical carrier of DNA, is packed inside of the nucleus in
complex ways we don't fully understand. We know that chromosomes tend to cluster
in characteristic locations, that certain regions of DNA tend to interact with
certain others (e.g. see [topologically associating domains][Beagan2020]), and
that certain interactions are crucial for biological function (e.g. see 
[promoter-enhancer interactions][Mora2015]). However, when these interactions
are studied, if a differential component is included, it is typically a
comparison between conditions or cell types. The goal of `bach` to elucidate the
**effect of genetic variation on chromatin interactions**.

#### Inputs

As input, `bach` uses read pairs from [Hi-C][Belton2012], a molecular biology
technique to detect chromatin segments in physical proximity. These read pairs
are assumed to come from **sample-specific [BAM][SAMspec] files** which have
already been sorted and indexed, for example by [samtools][samtools]. All BAM
files should be in a single directory (`--bam-dir`) and have a shared extension
(`--bam-ext`), such that filenames are formatted as `<dir>/<sample>.<ext>`.

Additionally, `bach` uses a [VCF][VCFspec] (`--vcf`) to know which variants to
test. The **sample names in the VCF must correspond to the BAM file names**.
The VCF may have extra samples, but it must have genotypes at least for all
samples with BAMs. Variants will be dropped from the VCF if they:
* do not have exactly two alleles,
* have alleles of lengths other than 1 (i.e. non-SNPs)
* failed filtration

#### Outputs

The allele-biased segments (see [Algorithm][#algorithm]) are written to a CSV
(`--output`) with these columns:
* `CHROM`: chromosome of SNP; VCF first column
* `SNP_POS`: 1-indexed bp position of SNP; VCF second column
* `SNP_ID`: ID of SNP; VCF third column
* `BIAS_START_MBP`: 0-indexed Mbp starting position of window with bias
* `BIAS_END_MBP`: 0-indexed Mbp ending position of window with bias
* `BIAS`: direction of bias (`ref`/`alt`)
* `P_VALUE`: p-value, from two-tailed binomial test, for this many samples
sharing the same bias

Additionally, a summary will be printed:
```
{A} found in test/bams/<sample>.pos.sorted.bam
{B} variants read from VCF
{C} non-filtered biallelic SNPs with >={D} heterozygotes tested
{E} other variants dropped
{F} biased segments written to CSV
```

`{A}` is the number of samples used in testing. `{C}` and `{E}` will add up to
`{B}`, which is itself the total number of variants within the VCF. `{D}` is the
minimum number of heterozyotes required to test a SNP for bias, `{A}` minus
`--max-drop`. `{E}` is the total number of variants dropped for any reason,
including both variants which are not non-filtered biallelic SNPs, as well as
SNPs which had too few heterozygotes to test. `{F}` is the total number of
biased segments found across all variants. 

> [!TIP]
> To also print a progress summary after every X variants processed, use
> `--print-progress X`. This is useful to track/estimate time for large jobs.

#### Algorithm

[Documentation under construction!]

[Beagan2020]: https://doi.org/10.1038/s41588-019-0561-1
[Belton2012]: https://doi.org/10.1016/j.ymeth.2012.05.001
[CondaInstall]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
[Chromatin]: https://www.genome.gov/genetics-glossary/Chromatin
[Deprecation]: https://docs.scipy.org/doc/scipy-1.11.1/reference/generated/scipy.stats.binom_test.html
[Mora2015]: https://doi.org/10.1093/bib/bbv097
[PGL]: https://pgl.soe.ucsc.edu/
[ProjectSummary]: https://docs.google.com/document/d/1V_6JWuLDxENUomvSkBamUB0kLm_riqedQ75qiWdRyVU/edit?usp=sharing
[pysamInstall]: https://niyunyun-pysam-fork.readthedocs.io/en/latest/installation.html#installationpytho
[PythonInstall]: https://www.python.org/downloads/
[SAMspec]: https://samtools.github.io/hts-specs/SAMv1.pdf
[samtools]: https://www.htslib.org/doc/samtools.html
[SciPyInstall]: https://scipy.org/install/
[VCFspec]: https://samtools.github.io/hts-specs/VCFv4.2.pdf
