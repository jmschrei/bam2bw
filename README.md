## bam2bw

[![Downloads](https://static.pepy.tech/badge/bam2bw)](https://pepy.tech/project/bam2bw)
[![Tests](https://github.com/jmschrei/bam2bw/actions/workflows/test.yml/badge.svg)](https://github.com/jmschrei/bam2bw/actions/workflows/test.yml)

A command-line tool for converting SAM/BAM files of reads, or .tsv/tsv.gz files of fragments, into either stranded or unstraded basepair resolution bigWig files. By default, only the 5' end of reads are mapped (not the full span of the read) and these bigWig file(s) contain the integer count of reads mapping to each basepair. Optionally, both the 3' and 5' of the entry can be mapped if they correspond to fragments, such as from ATAC-seq experiments. As a convenience, the starts and ends can be shifted (e.g., to account for Tn5 bias), a scaling factor can be used to multiply the mapped counts at each basepair, and read depth normalization can be applied to make the sum across the bigWigs be equal to 1. When a scaling factor and read depth normalization are used together, the sum across the two bigWigs is equal to the scaling factor.

`bam2bw` does not produce any intermediary files and can even stream SAM/BAM files remotely (but not .tsv/.tsv.gz). This means that you can go directly from finding a SAM/BAM file somewhere on the internet to the bigWig files used to train ML programs without several time-consuming steps. v0.4.0 allows parallel processing of files, even if they are remote, reducing the time needed to process inputs to just the time needed to process the biggest one.

```
usage: bam2bw [-h] -s SIZES [-u] [-f | -3p] [-ps POS_SHIFT] [-ns NEG_SHIFT] [-mp] [--rna5 {read1,read2}] [--opposite_strand] [-sf SCALE_FACTOR] [-r] [-p PARALLEL] -n NAME [-z ZOOMS] [-v] filename [filename ...]

This tool will convert BAM files to bigwig files without an intermediate.

positional arguments:
  filename              The SAM/BAM or tsv/tsv.gz file to be processed.

options:
  -h, --help            show this help message and exit
  -s SIZES, --sizes SIZES
                        A chrom_sizes, .fai, or FASTA file. Only the first two
                        columns of a chrom_sizes/.fai file are read. A
                        compressed FASTA must be BGZF, not gzip.
  -u, --unstranded      Have only one, unstranded, output.
  -f, --fragments       The data is fragments and so both ends should be recorded.
  -3p, --three_prime    Record the 3' end of each read instead of the 5' end.
  -ps POS_SHIFT, --pos_shift POS_SHIFT
                        A shift to apply to positive strand reads.
  -ns NEG_SHIFT, --neg_shift NEG_SHIFT
                        A shift to apply to negative strand reads.
  -mp, --mate_pairs     Treat paired-end BAM/SAM reads as a single RNA/fragment tag instead of
                        counting each mate independently (see --rna5/--opposite_strand).
  --rna5 {read1,read2}  Which mate carries the 5' end of the RNA/fragment. Only used with
                        --mate_pairs. Default: read1.
  --opposite_strand     Report the strand of the mate opposite the one chosen by --rna5.
                        Only used with --mate_pairs.
  -sf SCALE_FACTOR, --scale_factor SCALE_FACTOR
                        A scaling factor to multiply each position by.
  -r, --read_depth      Whether to divide through by total (pre-scaled) read depth.
  -p PARALLEL, --parallel PARALLEL
                        The number of jobs to use, max of one per input file.
  -n NAME, --name NAME
  -z ZOOMS, --zooms ZOOMS
                        The number of zooms to store in the bigwig.
  -v, --verbose
```

### Installation

`pip install bam2bw`

### Timings

These timings involve the processing of https://www.encodeproject.org/files/ENCFF638WXQ/ which has slightly over 70M reads. Local means applied to a file that was already downloaded, and remote means including the downloading time.

```
bam2bw (local): 2m10s
bam2bw (remote): 4m50s
existing pipeline (local): 18m5s
```

On my compute server, I can usually get over 1.5M records/second when reading BAM files locally and have been able to stream three BAMs from the ENCODE Portal at ~800k records/second, processing the almost 945M records in just over 8 minutes (ENCFF337YBN.bam, ENCFF981FXV.bam, ENCFF144GBU.bam). Admittedly, the first setting will be influenced by the quality of your hard drive, and the second setting by the quality of your internet connection (and whether the speed is throttled).

### Usage

(1) On a local file:

`bam2bw my.bam -s hg38.chrom.sizes -n test-run -v`

(2) On several local files:

`bam2bw my1.bam my2.bam my3.bam -s hg38.chrom.sizes -n test-run -v`

(3) On a remote file:

`bam2bw https://path/to/my.bam -s hg38.chrom.sizes -n test-run -v`

(4) On several remote files:

`bam2bw https://path/to/my1.bam https://path/to/my2.bam https://path/to/my3.bam -s hg38.chrom.sizes -n test-run -v`

Each will return two bigWig files: `test-run.+.bw` and `test-run.-.bw`. When multiple files are passed in their reads are concatenated without the need to produce an intermediary file of concatenated reads.

(5) When wanting a single unstranded bigWig:

`bam2bw my.bam -s hg38.chrom.sizes -n test-run -v -u`

(6) When wanting to map fragments (and get a single unstranded bigWig):

`bam2bw fragments.tsv.gz -s hg38.chrom.sizes -n test-run -v -f -u`

(7) With a FASTA instead of a .chrom.sizes:

`bam2bw my.bam -s hg38.fa -n test-run -v`

(8) When wanting to normalize by read depth such that the sum across both bigWigs is equal to 1.

`bam2bw my.bam -s hg38.chrom.sizes -n test-run -v -r`

(9) When wanting to normalize by read depth such that the sum across both bigWigs is equal to 1,000,000.

`bam2bw my.bam -s hg38.chrom.sizes -n test-run -v -r -sf 1000000`

(10) When wanting to record 3' ends instead of 5' ends:

`bam2bw my.bam -s hg38.chrom.sizes -n test-run -v -3p`

(11) When a BAM has paired-end reads that jointly represent a single RNA/fragment tag (e.g. PRO-seq/PRO-cap), rather than two independent events (e.g. the two Tn5 cut sites of an ATAC-seq fragment): use `--mate_pairs` so each pair contributes exactly one position instead of one from each mate. `--rna5` picks which mate carries the RNA's 5' end (the other mate's own 5' end is used as the RNA's 3' end), and `-3p`/`--opposite_strand` behave as before but are applied to the jointly-determined position/strand:

`bam2bw my.bam -s hg38.chrom.sizes -n test-run -v -mp --rna5 read2 -3p`

#### A note on paired-end BAMs

By default, `bam2bw` (like `bedtools genomecov`) counts every mapped alignment record independently, including both mates of a pair. For ATAC-seq/DNase-seq/ChIP-seq this is correct: each mate's end is its own real cut/fragment-boundary event. For assays where a fragment carries exactly one meaningful tag position determined jointly by both mates (PRO-seq, PRO-cap, and similar run-on/CAGE-style protocols), counting both mates independently will roughly double the signal and scatter it across the wrong positions (each mate maps to a different point in the fragment). Use `--mate_pairs` (see example 11) for that case instead.

#### Existing Pipeline

This tool is meant specifically to replace the following pipeline which produces several large intermediary files:

```bash
wget https://path/to/my.bam -O my.bam
samtools sort my.bam -o my.sorted.bam

bedtools genomecov -5 -bg -strand + -ibam my.sorted.bam | sort -k1,1 -k2,2n > my.+.bedGraph
bedtools genomecov -5 -bg -strand - -ibam my.sorted.bam | sort -k1,1 -k2,2n > my.-.bedGraph

bedGraphToBigWig my.+.bedGraph hg38.chrom.sizes my.+.bw
bedGraphToBigWig my.-.bedGraph hg38.chrom.sizes my.-.bw
```

### Testing

The test suite runs `bam2bw` end-to-end on small synthetic BAM, BED, and tsv
files built on the fly, and checks the values inside the bigWigs it writes.

```
pip install -e .[test]
pytest
```

or, with `uv`,

```
uv sync --extra test
uv run pytest
```

The same suite runs on GitHub Actions on every push and pull request to `main`,
across Python 3.10 through 3.13.

### Version Log

```
v0.5.0
======

  - Added -3p/--three_prime to record the 3' end of each read or interval instead of the 5' end.
  - Added -mp/--mate_pairs, --rna5, and --opposite_strand to jointly count paired-end
    reads as a single RNA/fragment tag (e.g. for PRO-seq/PRO-cap) instead of counting
    each mate independently.
v0.4.3
======

  - Packaging moved from setup.py to pyproject.toml. `pip install bam2bw` and `pip install -e .[test]` are unchanged.
  - Entries falling outside the chromosome sizes given to -s are now discarded explicitly, with a count reported, instead of being dropped by pyBigWig without notice.
  - -r now normalizes over the entries actually written rather than over every counted read, so a sizes file that disagrees with the input no longer yields a track summing to less than the scale factor.
  - -r with -v reports the read depth that was divided through.
  - A samtools .fai index, or any sizes file with more than two columns, can now be passed to -s; only the first two fields of each line are read.
  - Blank lines and # comments in a chrom_sizes file are skipped instead of raising.
  - Malformed BED/tsv lines and mapped reads with no CIGAR now report the file and the offending line or read rather than raising a bare unpacking error.
  - A gzip-compressed FASTA passed to -s now says that BGZF is required and names the file. Only BGZF has ever worked; the v0.4.2 note about .gz variants refers to BGZF.

v0.4.2
======

  - Recognize additional FASTA extensions (.fasta, .fna, .fas) and their .gz variants when passed to -s.
  - Fixed tqdm progress bars clashing when processing input files in parallel.
  - Added biopython as an explicit dependency.

v0.4.1
======

  - Oops removed print statement.

v0.4.0
======

  - Added in a -p/--parallel option to read from input files in parallel using joblib. Max 1 process per input file
  - Moved the opening of bigWig objects until AFTER the input files are completely read to avoid deleting everything immediately if you accidentally re-run the command
  - Allow reading from .bed and .bed.gz files.
  - Create one progress bar per input file that can update in parallel with the description being the file name.


v0.3.2
======

  - Added -sf which is a scale factor that multiplies each position.
  - Added -r which performs read-depth normalization, dividing each position by the (pre-scaled) sum.
  - You can use -sf and -r together to make the bigWigs sum to your desired value.


v0.3.1
======

  - Fixed a minor bug

v0.3.0
======

  - A FASTA can be passed in to -s instead of a chrom_sizes file, and the chromosomes and their sizes automatically extracted
  - Entries not mapping to a chromosome in the chrom_sizes/FASTA file will be ignored and a warning will be raised if -v is set
  - .tsv and .tsv.gz files can processed now using the same arguments.
  - The -f argument has been added which treats entires as fragments where both the 3' and 5' ends should be added instead of just the 5' one

```
