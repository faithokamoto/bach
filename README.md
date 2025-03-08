# `bach`: Biased Allele Contacts in Hi-C

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

Alternatively, run [conda][CondaInstall] from inside of the `bach` directory:
```
conda env create --file environment.yml
conda activate bach
```

## Examples/tests

To run `bach` on the provided test data, use 
```
./bach --max-opposite 0 --max-neutral 0 --max-drop 1 \
    --window-width 15 --window-step 5 \
    -v test/snps.vcf -d test/bams -o output.csv
```

This commands correspond to: "Using BAMs in `test/bams` and VCF `test/snps.vcf`,
output biased segments to `output.csv`. At most one sample per SNP may be
dropped/ignored due to a missing or homozygous genotype. No samples within a 
biased window may have opposite or neutral bias. Use windows of 15Mbp with a
step of 5Mbp." More algorithm details are in [Algorithm](#algorithm).

The output files should be the same as `test/output.csv`

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
            [-N 0.*] [-D #] [-w Mbp] [-s Mbp]

Find allele-biased segments. Compare mates of reads with SNP alleles.

options:
  -h, --help            show this help message and exit
  -p #, --progress #    Progress update frequency (# variants processed)

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
                        Minimum window width, in Mbp, to scan with
  -s Mbp, --window-step Mbp
                        Window step size, in Mbp
```

### Background

[Chromatin], the physical carrier of DNA, is packed in the nucleus. It is known
that chromosomes cluster in certain locations, that certain regions interact
with certain others (e.g. see [topologically associating domains][Beagan2020]),
and that certain interactions are crucial for biological function (e.g. see 
[promoter-enhancer interactions][Mora2015]). However, when these interactions
are studied, if a differential component is included, it is typically a
comparison between conditions or cell types. The goal of `bach` to elucidate the
**effect of genetic variation on chromatin interactions**.

### Inputs

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

### Outputs

The allele-biased segments (see [Algorithm](#algorithm)) are written to a CSV
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
> To also print a progress summary after every `X` variants processed, use
> `--print-progress X`. This is useful to track/estimate time for large jobs.

### Problem formulation

The goal of `bach` is to find **allele-biased chromatin contacts**. These are
defined as regions of a chromosome which preferentially contact one allele over
another, as measured by Hi-C read mates. That is, given:
* *n* samples
* *O* ∈ [0, 0.5) maximum share of heterozygous samples with opposite bias
* *N* ∈ [0, 1) maximum share of heterozygous samples with neutral bias
* *D* ∈ [0, *n*) maximum number of samples which may be dropped
* *G*<sub>*i*</sub> sample genotypes for ∀ *i* ∈ [1, *n*]
    * *k* ⊆ [1, *n*] where *i* ∈ *k* if *G*<sub>*i*</sub> is heterozygous
    * **|*k*| ≥ *n* − *D* is assumed**
* *M*<sub>*a*,*i*</sub> read mate sets for ∀ *i* ∈ [1, *n*] and for ∀ *a* ∈
[ref, alt] for the positions of all mates of sample *i*'s reads with allele *a*
* *p* window location precision
* *w* minimum window width

Find all windows [*x*, *y*] (for *x* and *y* multiples of *p* and *y* − *x* ≥
*w*) on the variant's chromosome such that for some allele *a*:
* given mate counts by allele:
    * *C*<sub>*A*</sub> = size( *M*<sub>*a*,*i*</sub> ∈ [*x*, *y*] ) for sample *i*
    * *C*<sub>*O*</sub> = size( *M*<sub>*~a*,*i*</sub> ∈ [*x*, *y*] ) for sample *i*
* and thus sample bias counts:
    * *B*<sub>*A*</sub> = Σ<sub>*i* ∈ *k*</sub> ( *C*<sub>*A*</sub> > *C*<sub>*O*</sub> )
    * *B*<sub>*O*</sub> = Σ<sub>*i* ∈ *k*</sub> ( *C*<sub>*A*</sub> < *C*<sub>*O*</sub> )
    * *B*<sub>*N*</sub> = Σ<sub>*i* ∈ *k*</sub> ( *C*<sub>*A*</sub> = *C*<sub>*O*</sub> )
* the window is biased as so:
    * *B*<sub>*A*</sub> > *B*<sub>*O*</sub>
    * *B*<sub>*O*</sub> ≤ |*k*| × *O*
    * *B*<sub>*N*</sub> ≤ |*k*| × *N*

As a concrete example, take *n*=11, *O*=0.1, *N*=0.2, and *D*=2. At most 2/11
samples may be dropped due to a missing or homozygous genotype. Assume one was
dropped, leaving ten heterozygous samples for testing. At most `10 × 0.1 = 1`
may have a bias opposite of the majority, and at most `10 × 0.2 = 2` may have
neutral/no bias. That is, in the worst case to call a biased window, at least
`10 - 2 - 1 = 7` samples must have the same bias (e.g., favoring the reference
allele), one may have the opposite bias (e.g., favoring the alternate), and two
may have no bias.

To determine if a specific one of those samples is biased in a specific window,
simply count mates of reads with a given allele. Say that there are four reads
pairs where one mate overlaps the variant and one mate is within the window. If
two of those mates have allele *a*, while the other two have the other allele,
then this window has no bias. On the other hand, if all of those mates have
allele *a*, that indicates that haplotypes with allele *a* are more likely to
contact the chromosome within this window, and thus it is biased towards *a*.

> [!NOTE]
> While the above formulation is agnostic towards the type of variant, `bach`
> only supports single nucleotide polymorphism ([SNP][SNP]) variants as of now.
> Also, while it is possible to formulate this problem more generally to support
> variants with more than two alleles, `bach` only supports biallelic SNPs.
> Finally, note that `bach` ignores contacts to other chromosomes.

### Algorithm

The general idea of `bach`'s algorithm is to do the following *for each SNP*:
1. Subset samples to only those which are **heterozygous**. If more than
`--max-drop` samples have homozygous or missing genotypes, then skip this SNP.
2. Look up all reads overlapping the SNP. Save the locations of their **mates**
in two lists, separated by the allele of the mate overlapping the SNP.
3. Slide a **window** (width `--window-width`) across the SNP's chromosome.
4. For each window, determine *for each sample* whether mates of one allele
are more abundant than the other.
    * Windows which **overlap the SNP** itself are skipped.
    * Window move in steps of `--move-step` Mbp, starting from 0, until the end
    of the chromosome.
    * This portion of the algorithm uses pointers to indices in the sorted lists
    of mate positions, and thus will only traverse each list once per variant.
5. Determine whether this window is biased by **counting** ref-biased samples,
alt-biased samples, and neutral (identical mate count for both alleles) samples.
    * The overall bias of the window is either `ref` or `alt` based on which
    category has more biased samples.
    * If more than `--max-opposite` share of the samples have the opposite bias
    (e.g. `alt` when more samples are `ref` biased), drop this window.
    * If more than `--max-neutral` share of the samples have a neutral bias,
    drop this window.

> [!WARNING]
> Unlike `--max-drop`, an integer absolute number of samples which may be
> dropped, both `--max-opposite` and `--max-neutral` are **decimal shares**
> of the total number of heterozygotes. If there are 20 total samples, but
> two of them were dropped due to a homozygous genotype, then the number
> which may have a neutral bias is `(20 - 2) × --max-neutral`.
6. Assuming this window wasn't dropped:
    1. Attempt to extend it by `p`. That is, check a mini-window of length `p`
    starting from the end of the current window. If the extension has the same
    bias, add it to the window. Repeat until extension is no longer possible.
    2. Calculate a **two-tailed [binomial][BinomTest] p-value**, with `n` as the
    total number of samples with either `ref` or `alt` bias (i.e. ignoring
    neutral samples), `k` as the number of samples with majority bias, and `p`
    as `0.5` (i.e. with the null hypothesis that either bias is equally likely).

[BinomTest]: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binomtest.html
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
[SNP]: https://www.nature.com/scitable/definition/snp-295/
[VCFspec]: https://samtools.github.io/hts-specs/VCFv4.2.pdf
