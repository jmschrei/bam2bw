## bam2bw

[![Downloads](https://static.pepy.tech/badge/bam2bw)](https://pepy.tech/project/bam2bw)

A command-line tool for converting SAM/BAM files of reads, or .tsv/tsv.gz files of fragments, into either stranded or unstraded basepair resolution bigWig files. By default, only the 5' end of reads are mapped (not the full span of the read) and these bigWig file(s) contain the integer count of reads mapping to each basepair. Optionally, both the 3' and 5' of the entry can be mapped if they correspond to fragments, such as from ATAC-seq experiments. As a convenience, the starts and ends can be shifted (e.g., to account for Tn5 bias), a scaling factor can be used to multiply the mapped counts at each basepair, and read depth normalization can be applied to make the sum across the bigWigs be equal to 1. When a scaling factor and read depth normalization are used together, the sum across the two bigWigs is equal to the scaling factor.

`bam2bw` does not produce any intermediary files and can even stream SAM/BAM files remotely (but not .tsv/.tsv.gz). This means that you can go directly from finding a SAM/BAM file somewhere on the internet to the bigWig files used to train ML programs without several time-consuming steps.

```
usage: bam2bw [-h] -s SIZES [-u] [-f] [-ps POS_SHIFT] [-ns NEG_SHIFT] [-sf SCALE_FACTOR] [-r] -n NAME [-z ZOOMS] [-v] filename [filename ...]

This tool will convert BAM files to bigwig files without an intermediate.

positional arguments:
  filename              The SAM/BAM or tsv/tsv.gz file to be processed.

options:
  -h, --help            show this help message and exit
  -s SIZES, --sizes SIZES
                        A chrom_sizes or FASTA file.
  -u, --unstranded      Have only one, unstranded, output.
  -f, --fragments       The data is fragments and so both ends should be recorded.
  -ps POS_SHIFT, --pos_shift POS_SHIFT
                        A shift to apply to positive strand reads.
  -ns NEG_SHIFT, --neg_shift NEG_SHIFT
                        A shift to apply to negative strand reads.
  -sf SCALE_FACTOR, --scale_factor SCALE_FACTOR
                        A scaling factor to multiply each position by.
  -r, --read_depth      Whether to divide through by total (pre-scaled) read depth.
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

### Version Log

```
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
