## bam2bw

A command-line tool for converting SAM/BAM files into either stranded or unstranded basepair resolution bigWig files. Specifically, only the 5' end of reads are mapped (not the full span of the read) and these bigWig file(s) contain the integer count of reads mapping to each basepair. `bam2bw` does not produce any intermediary files and can even stream SAM/BAM files remotely. This means that you can go directly from finding a SAM/BAM file somewhere on the internet to the BigWig files used to train ML programs without several time-consuming steps.

```
usage: bam2bw [-h] -s SIZES [-u] -n NAME [-ps POS_SHIFT] [-ns NEG_SHIFT] [-z ZOOMS] [-v] filename [filename ...]

This tool will convert BAM files to bigwig files without an intermediate.

positional arguments:
  filename              The SAM/BAM file to be processed.

options:
  -h, --help            show this help message and exit
  -s SIZES, --sizes SIZES
                        A chromosome sizes file.
  -u, --unstranded      Have only one, unstranded, output.
  -n NAME, --name NAME
  -ps POS_SHIFT, --pos_shift POS_SHIFT
                        A shift to apply to positive strand reads.
  -ns NEG_SHIFT, --neg_shift NEG_SHIFT
                        A shift to apply to negative strand reads.
  -z ZOOMS, --zooms ZOOMS
                        The number of zooms to store in the bigwig.
  -v, --verbose
```

#### Installation

`pip install bam2bw`

#### Timings

These timings involve the processing of https://www.encodeproject.org/files/ENCFF638WXQ/ which has slightly over 70M reads. Local means applied to a file that was already downloaded, and remote means including the downloading time.

```
bam2bw (local): 2m10s
bam2bw (remote): 4m50s
existing pipeline (local): 18m5s
```

#### Usage

On a local file:

`bam2bw my.bam -s hg38.chrom.sizes -n test-run -v`

On several local files:

`bam2bw my1.bam my2.bam my3.bam -s hg38.chrom.sizes -n test-run -v`

On a remote file:

`bam2bw https://path/to/my.bam -s hg38.chrom.sizes -n test-run -v`

On several remote files:

`bam2bw https://path/to/my1.bam https://path/to/my2.bam https://path/to/my3.bam -s hg38.chrom.sizes -n test-run -v`

Each will return two BigWig files: `test-run.+.bw` and `test-run.-.bw`. When multiple files are passed in their reads are concatenated without the need to produce an intermediary file of concatenated reads.

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

