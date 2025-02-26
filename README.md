# bach
Biased Allele Contacts in Hi-C

BME 230A class project, heavily based on a [rotation project][ProjectSummary]
performed in [Dr. Richard "Ed" Green's lab][PGL] working with Dr. Merly Escalona.

## Installation

Clone `bach` from GitHub:

```
git clone https://github.com/faithokamoto/bach
cd bach # Enter directory
chmod +x ./bach # Make executable
```

### Dependencies

`bach` has the following dependencies, documented in `environment.yml`:
* Python 3 (tested on version `3.12.8`) [Installation instructions][PythonInstall]
* `pysam` (tested on version `0.23.0`) [Installation instructions][pysamInstall]
* `SciPy` (tested on version `1.15.2`) [Installation instructions][ScipyInstall]

> [!NOTE]
> Of these dependencies, the one most likely to have an issue is `SciPy`.
> Versions before 1.10.0 will lack `scipy.stats.binomtest`, instead using
> `scipy.stats.binom_test`, which has since been [deprecated][Deprecation].

Alternatively, if you have [conda][CondaInstall], simply run this from inside of
the `bach` directory:
```
conda env create --file environment.yml
```

Test by:
```
./bach --help
```

## Test data

To run `bach` on the provided test data, use 
```
./bach --max-opposite 0 --max-neutral 0 --max-drop 1 test/snps.vcf test/bams test/out.csv
```

To generate new test data, use 
```
python scripts/simulate_reads.py --depth 30 \
    --biased NC_000020.11:1182052 --biased NC_000020.11:1182049 \
    test/sequence.fasta test/snps.vcf test/bams
```

## Options

[Documentation in progress!]


[CondaInstall]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html
[Deprecation]: https://docs.scipy.org/doc/scipy-1.11.1/reference/generated/scipy.stats.binom_test.html
[PGL]: https://pgl.soe.ucsc.edu/
[ProjectSummary]: https://docs.google.com/document/d/1V_6JWuLDxENUomvSkBamUB0kLm_riqedQ75qiWdRyVU/edit?usp=sharing
[pysamInstall]: https://niyunyun-pysam-fork.readthedocs.io/en/latest/installation.html#installationpytho
[PythonInstall]: https://www.python.org/downloads/
[SciPyInstall]: https://scipy.org/install/
