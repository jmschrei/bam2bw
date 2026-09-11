# conftest.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com>

import os
import subprocess
import sys

import pytest

from .bigwig import read_bigwig
from .synthetic import write_bam
from .synthetic import write_chrom_sizes
from .synthetic import write_fasta
from .synthetic import write_intervals


BAM2BW = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
	"bam2bw")


# The one read layout that almost every test in the suite runs against. It is
# deliberately small enough to compute expected bigWig values by hand, while
# still covering each branch the main loop can take:
#
#    - both strands, on more than one chromosome,
#    - positions hit by more than one read, so that counts have to accumulate
#      rather than overwrite,
#    - a chromosome in the sizes file that no read touches (chr3),
#    - a chromosome the reads touch that is not in the sizes file (chrUN),
#    - an unmapped read.
#
# The 5' end of a forward read is its reference start and its 3' end is
# reference_end - 1; for a reverse read the two are swapped. Working those out
# for the reads below gives, with no flags at all:
#
#    plus    chr1 {100: 2, 200: 1}    chr2 {50: 1}
#    minus   chr1 {339: 2, 419: 1}    chr2 {89: 1}
#
# for a total of eight counted reads, which is a power of two so that the
# read-depth normalized values stay exactly representable as float32 and can
# still be compared to four decimal places.

CHROM_SIZES = [("chr1", 1000), ("chr2", 500), ("chr3", 200)]
BAM_CHROM_SIZES = CHROM_SIZES + [("chrUN", 100)]

READS = [
	("chr1", 100, 50, False),
	("chr1", 100, 50, False),
	("chr1", 200, 30, False),
	("chr1", 300, 40, True),
	("chr1", 300, 40, True),
	("chr1", 400, 20, True),
	("chr2", 50, 25, False),
	("chr2", 80, 10, True),
	("chrUN", 10, 10, False)
]

# The unmapped read is placed on chr1, the way a read whose mate aligned is
# placed in a real BAM. It must be dropped by the unmapped check rather than
# by the missing-chromosome check, so that removing either one shows up.

UNMAPPED = [("chr1", 500, 20)]

# The interval files carry no strand, so `bam2bw` treats every entry as being
# on the forward strand. The layout mirrors the forward-strand reads above so
# that the two input paths can be compared against each other.

INTERVALS = [
	("chr1", 100, 150),
	("chr1", 100, 150),
	("chr1", 200, 230),
	("chr2", 50, 75),
	("chrUN", 10, 20)
]

FASTA_EXTENSIONS = ".fa", ".fasta", ".fna", ".fas"


## Input files, built once for the whole session


@pytest.fixture(scope="session")
def data_dir(tmp_path_factory):
	return tmp_path_factory.mktemp("data")


@pytest.fixture(scope="session")
def sizes(data_dir):
	return write_chrom_sizes(data_dir / "test.chrom.sizes", CHROM_SIZES)


@pytest.fixture(scope="session")
def fastas(data_dir):
	"""A FASTA of the same chromosomes under each accepted extension."""

	return {extension: write_fasta(data_dir / ("test" + extension), CHROM_SIZES)
		for extension in FASTA_EXTENSIONS}


@pytest.fixture(scope="session")
def bam(data_dir):
	return write_bam(data_dir / "test.bam", BAM_CHROM_SIZES, READS,
		unmapped=UNMAPPED)


@pytest.fixture(scope="session")
def bed(data_dir):
	return write_intervals(data_dir / "test.bed", INTERVALS)


@pytest.fixture(scope="session")
def bed_gz(data_dir):
	return write_intervals(data_dir / "test.bed.gz", INTERVALS)


@pytest.fixture(scope="session")
def tsv(data_dir):
	return write_intervals(data_dir / "test.tsv", INTERVALS)


@pytest.fixture(scope="session")
def tsv_gz(data_dir):
	return write_intervals(data_dir / "test.tsv.gz", INTERVALS)


## Running the tool


@pytest.fixture
def run(tmp_path):
	"""Run `bam2bw` in a scratch directory and hand back the finished process.

	The `-n` prefix is supplied automatically and points into this test's own
	`tmp_path`, so tests never name output files and never collide. Nothing is
	asserted about the return code here -- use this fixture directly only when
	the failure itself is what is being tested.
	"""

	def _run(*args, name="out"):
		command = [sys.executable, BAM2BW]
		command += [str(arg) for arg in args]
		command += ["-n", str(tmp_path / name)]

		return subprocess.run(command, capture_output=True, text=True)

	return _run


@pytest.fixture
def stranded(tmp_path, run):
	"""Run `bam2bw` and read back the plus and minus bigWigs it wrote."""

	def _stranded(*args, name="out"):
		process = run(*args, name=name)

		assert process.returncode == 0, process.stderr

		return (read_bigwig(tmp_path / (name + ".+.bw")),
			read_bigwig(tmp_path / (name + ".-.bw")))

	return _stranded


@pytest.fixture
def unstranded(tmp_path, run):
	"""Run `bam2bw -u` and read back the single bigWig it wrote."""

	def _unstranded(*args, name="out"):
		process = run(*args, "-u", name=name)

		assert process.returncode == 0, process.stderr

		return read_bigwig(tmp_path / (name + ".bw"))

	return _unstranded
