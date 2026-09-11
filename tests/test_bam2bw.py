# test_bam2bw.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com>

import subprocess
import sys

import pysam
import pytest

from numpy.testing import assert_array_almost_equal

from .bigwig import entries
from .bigwig import read_bigwig
from .bigwig import total

from .conftest import BAM2BW
from .conftest import CHROM_SIZES
from .conftest import FASTA_EXTENSIONS

from .synthetic import write_chrom_sizes
from .synthetic import write_intervals


# `bam2bw` is a script rather than an importable module, so every test here
# runs the real command line in a subprocess and reads the bigWigs back. That
# is also the contract users actually depend on: the flags, the two-or-one
# output files, and the numbers inside them.
#
# The read layout and the expected values it implies are documented in
# conftest.py. Expected positions and counts are hardcoded below rather than
# recomputed, so that a change in the recording logic shows up as a failing
# test instead of as two matching re-derivations of the same mistake.


## Input validation


@pytest.mark.parametrize("extension", [".txt", ".fastq", ".bw", ".bigWig",
	".bam.bai", ".fa", ""])
def test_rejects_unknown_extension(run, sizes, tmp_path, extension):
	path = tmp_path / ("reads" + extension)
	path.write_text("")

	process = run(path, "-s", sizes)

	assert process.returncode != 0
	assert "Filenames must end in one of" in process.stderr


def test_rejects_unknown_extension_among_valid_ones(run, sizes, bam, tmp_path):
	path = tmp_path / "reads.txt"
	path.write_text("")

	process = run(bam, path, "-s", sizes)

	assert process.returncode != 0
	assert "Filenames must end in one of" in process.stderr


def test_requires_sizes(run, bam):
	process = run(bam)

	assert process.returncode == 2
	assert "-s/--sizes" in process.stderr


def test_requires_name(bam, sizes):
	# The -n prefix is supplied by the `run` fixture, so this one test has to
	# build the command itself in order to leave it out.
	process = subprocess.run([sys.executable, BAM2BW, str(bam), "-s",
		str(sizes)], capture_output=True, text=True)

	assert process.returncode == 2
	assert "-n/--name" in process.stderr


def test_requires_an_input_file(run, sizes):
	process = run("-s", sizes)

	assert process.returncode == 2


def test_missing_input_file(run, sizes, tmp_path):
	process = run(tmp_path / "does_not_exist.bam", "-s", sizes)

	assert process.returncode != 0


def test_missing_sizes_file(run, bam, tmp_path):
	process = run(bam, "-s", tmp_path / "does_not_exist.chrom.sizes")

	assert process.returncode != 0


## Chromosome sizes and FASTA input


def test_chrom_sizes_sets_the_header(stranded, bam, sizes):
	pos, neg = stranded(bam, "-s", sizes)

	assert list(pos.keys()) == ["chr1", "chr2", "chr3"]
	assert list(neg.keys()) == ["chr1", "chr2", "chr3"]


def test_chrom_not_in_sizes_is_absent_from_output(stranded, bam, sizes):
	pos, neg = stranded(bam, "-s", sizes)

	assert "chrUN" not in pos
	assert "chrUN" not in neg


@pytest.mark.parametrize("extension", FASTA_EXTENSIONS)
def test_fasta_is_accepted(stranded, bam, fastas, extension):
	pos, neg = stranded(bam, "-s", fastas[extension])

	assert list(pos.keys()) == ["chr1", "chr2", "chr3"]

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 200])
	assert_array_almost_equal(counts, [2, 1], 4)


@pytest.mark.parametrize("extension", FASTA_EXTENSIONS)
def test_fasta_matches_chrom_sizes(stranded, bam, sizes, fastas, extension):
	from_sizes, _ = stranded(bam, "-s", sizes, name="sizes")
	from_fasta, _ = stranded(bam, "-s", fastas[extension], name="fasta")

	assert from_sizes == from_fasta


def test_chrom_sizes_order_is_preserved(stranded, bam, tmp_path):
	reordered = write_chrom_sizes(tmp_path / "reordered.chrom.sizes",
		[("chr3", 200), ("chr1", 1000), ("chr2", 500)])

	pos, _ = stranded(bam, "-s", reordered)

	assert list(pos.keys()) == ["chr3", "chr1", "chr2"]


## Output file naming


def test_stranded_writes_two_files(run, bam, sizes, tmp_path):
	process = run(bam, "-s", sizes)

	assert process.returncode == 0
	assert (tmp_path / "out.+.bw").exists()
	assert (tmp_path / "out.-.bw").exists()
	assert not (tmp_path / "out.bw").exists()


def test_unstranded_writes_one_file(run, bam, sizes, tmp_path):
	process = run(bam, "-s", sizes, "-u")

	assert process.returncode == 0
	assert (tmp_path / "out.bw").exists()
	assert not (tmp_path / "out.+.bw").exists()
	assert not (tmp_path / "out.-.bw").exists()


## BAM: 5' ends, the default behaviour


def test_bam_forward_reads(stranded, bam, sizes):
	pos, _ = stranded(bam, "-s", sizes)

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 200])
	assert_array_almost_equal(counts, [2, 1], 4)

	positions, counts = entries(pos, "chr2")
	assert_array_almost_equal(positions, [50])
	assert_array_almost_equal(counts, [1], 4)


def test_bam_reverse_reads(stranded, bam, sizes):
	_, neg = stranded(bam, "-s", sizes)

	positions, counts = entries(neg, "chr1")
	assert_array_almost_equal(positions, [339, 419])
	assert_array_almost_equal(counts, [2, 1], 4)

	positions, counts = entries(neg, "chr2")
	assert_array_almost_equal(positions, [89])
	assert_array_almost_equal(counts, [1], 4)


def test_bam_chromosome_without_reads_is_empty(stranded, bam, sizes):
	pos, neg = stranded(bam, "-s", sizes)

	assert pos["chr3"] == {}
	assert neg["chr3"] == {}


def test_bam_drops_unmapped_and_filtered_reads(stranded, bam, sizes):
	pos, neg = stranded(bam, "-s", sizes)

	# Nine reads plus one unmapped read go in; the unmapped read and the one
	# on chrUN are dropped, leaving eight.
	assert total(pos) + total(neg) == 8


## -f/--fragments


def test_bam_fragments_forward(stranded, bam, sizes):
	pos, _ = stranded(bam, "-s", sizes, "-f")

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 149, 200, 229])
	assert_array_almost_equal(counts, [2, 2, 1, 1], 4)

	positions, counts = entries(pos, "chr2")
	assert_array_almost_equal(positions, [50, 74])
	assert_array_almost_equal(counts, [1, 1], 4)


def test_bam_fragments_reverse(stranded, bam, sizes):
	_, neg = stranded(bam, "-s", sizes, "-f")

	positions, counts = entries(neg, "chr1")
	assert_array_almost_equal(positions, [300, 339, 400, 419])
	assert_array_almost_equal(counts, [2, 2, 1, 1], 4)

	positions, counts = entries(neg, "chr2")
	assert_array_almost_equal(positions, [80, 89])
	assert_array_almost_equal(counts, [1, 1], 4)


def test_bam_fragments_double_the_signal(stranded, bam, sizes):
	pos, neg = stranded(bam, "-s", sizes)
	pos_f, neg_f = stranded(bam, "-s", sizes, "-f", name="frag")

	assert total(pos_f) + total(neg_f) == 2 * (total(pos) + total(neg))


## -u/--unstranded


def test_unstranded_merges_both_strands(unstranded, bam, sizes):
	values = unstranded(bam, "-s", sizes)

	positions, counts = entries(values, "chr1")
	assert_array_almost_equal(positions, [100, 200, 339, 419])
	assert_array_almost_equal(counts, [2, 1, 2, 1], 4)

	positions, counts = entries(values, "chr2")
	assert_array_almost_equal(positions, [50, 89])
	assert_array_almost_equal(counts, [1, 1], 4)


def test_unstranded_equals_sum_of_strands(stranded, unstranded, bam, sizes):
	pos, neg = stranded(bam, "-s", sizes, name="stranded")
	values = unstranded(bam, "-s", sizes, name="unstranded")

	for chrom, _ in CHROM_SIZES:
		merged = dict(pos[chrom])
		for position, count in neg[chrom].items():
			merged[position] = merged.get(position, 0) + count

		assert values[chrom] == merged


def test_unstranded_fragments(unstranded, bam, sizes):
	values = unstranded(bam, "-s", sizes, "-f")

	positions, counts = entries(values, "chr1")
	assert_array_almost_equal(positions, [100, 149, 200, 229, 300, 339, 400,
		419])
	assert_array_almost_equal(counts, [2, 2, 1, 1, 2, 2, 1, 1], 4)

	positions, counts = entries(values, "chr2")
	assert_array_almost_equal(positions, [50, 74, 80, 89])
	assert_array_almost_equal(counts, [1, 1, 1, 1], 4)


## -ps/--pos_shift and -ns/--neg_shift


def test_zero_shifts_are_the_default(stranded, bam, sizes):
	pos, neg = stranded(bam, "-s", sizes, name="default")
	pos_z, neg_z = stranded(bam, "-s", sizes, "-ps", 0, "-ns", 0, name="zero")

	assert pos == pos_z
	assert neg == neg_z


def test_pos_shift_moves_only_the_recorded_start(stranded, bam, sizes):
	pos, neg = stranded(bam, "-s", sizes, "-ps", 4)

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [104, 204])
	assert_array_almost_equal(counts, [2, 1], 4)

	# Reverse reads record reference_end - 1, which -ps does not touch.
	positions, counts = entries(neg, "chr1")
	assert_array_almost_equal(positions, [339, 419])
	assert_array_almost_equal(counts, [2, 1], 4)


def test_neg_shift_moves_only_the_recorded_end(stranded, bam, sizes):
	pos, neg = stranded(bam, "-s", sizes, "-ns", 5)

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 200])
	assert_array_almost_equal(counts, [2, 1], 4)

	positions, counts = entries(neg, "chr1")
	assert_array_almost_equal(positions, [344, 424])
	assert_array_almost_equal(counts, [2, 1], 4)


def test_negative_shifts_are_accepted(stranded, bam, sizes):
	"""The Tn5 correction for ATAC-seq is +4 on one strand and -5 on the other.

	This is the single most common way `bam2bw` is invoked on ATAC-seq, and a
	negative value is easy to break at the argparse level without breaking
	anything else, so it gets its own test rather than only appearing inside a
	combination.
	"""

	pos, neg = stranded(bam, "-s", sizes, "-ps", 4, "-ns", -5)

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [104, 204])
	assert_array_almost_equal(counts, [2, 1], 4)

	positions, counts = entries(neg, "chr1")
	assert_array_almost_equal(positions, [334, 414])
	assert_array_almost_equal(counts, [2, 1], 4)

	positions, counts = entries(pos, "chr2")
	assert_array_almost_equal(positions, [54])
	assert_array_almost_equal(counts, [1], 4)

	positions, counts = entries(neg, "chr2")
	assert_array_almost_equal(positions, [84])
	assert_array_almost_equal(counts, [1], 4)


@pytest.mark.parametrize("flag", ["-ps", "-ns", "--pos_shift", "--neg_shift"])
def test_negative_shift_parses_on_every_spelling(run, bam, sizes, flag):
	process = run(bam, "-s", sizes, flag, -5)

	assert process.returncode == 0, process.stderr


def test_shifts_with_fragments(stranded, bam, sizes):
	pos, neg = stranded(bam, "-s", sizes, "-f", "-ps", 4, "-ns", -5)

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [104, 144, 204, 224])
	assert_array_almost_equal(counts, [2, 2, 1, 1], 4)

	positions, counts = entries(neg, "chr1")
	assert_array_almost_equal(positions, [304, 334, 404, 414])
	assert_array_almost_equal(counts, [2, 2, 1, 1], 4)

	# The chr2 reverse read is only 10bp long, so shifting both of its ends
	# inward collapses them onto the same position, where the counts add
	# rather than overwrite.
	positions, counts = entries(neg, "chr2")
	assert_array_almost_equal(positions, [84])
	assert_array_almost_equal(counts, [2], 4)


def test_shifts_preserve_total_signal(stranded, bam, sizes):
	pos, neg = stranded(bam, "-s", sizes, name="default")
	pos_s, neg_s = stranded(bam, "-s", sizes, "-ps", 4, "-ns", -5,
		name="shifted")

	assert total(pos_s) + total(neg_s) == total(pos) + total(neg)


## -sf/--scale_factor


def test_scale_factor_one_is_the_default(stranded, bam, sizes):
	pos, neg = stranded(bam, "-s", sizes, name="default")
	pos_s, neg_s = stranded(bam, "-s", sizes, "-sf", 1, name="scaled")

	assert pos == pos_s
	assert neg == neg_s


def test_scale_factor_multiplies_every_value(stranded, bam, sizes):
	pos, neg = stranded(bam, "-s", sizes, "-sf", 10)

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 200])
	assert_array_almost_equal(counts, [20, 10], 4)

	positions, counts = entries(neg, "chr1")
	assert_array_almost_equal(positions, [339, 419])
	assert_array_almost_equal(counts, [20, 10], 4)


def test_fractional_scale_factor(stranded, bam, sizes):
	pos, _ = stranded(bam, "-s", sizes, "-sf", 0.5)

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 200])
	assert_array_almost_equal(counts, [1, 0.5], 4)


def test_scale_factor_does_not_move_positions(stranded, bam, sizes):
	pos, _ = stranded(bam, "-s", sizes, name="default")
	pos_s, _ = stranded(bam, "-s", sizes, "-sf", 7.5, name="scaled")

	assert sorted(pos["chr1"]) == sorted(pos_s["chr1"])


## -r/--read_depth


def test_read_depth_normalizes_to_one(stranded, bam, sizes):
	pos, neg = stranded(bam, "-s", sizes, "-r")

	assert_array_almost_equal(total(pos) + total(neg), 1.0, 4)

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 200])
	assert_array_almost_equal(counts, [0.25, 0.125], 4)


def test_read_depth_counts_both_strands_once(unstranded, stranded, bam, sizes):
	"""Read depth is summed over both strands when stranded and over the one
	merged dictionary when not, so the two must agree rather than the
	unstranded run halving or doubling."""

	pos, neg = stranded(bam, "-s", sizes, "-r", name="stranded")
	values = unstranded(bam, "-s", sizes, "-r", name="unstranded")

	assert_array_almost_equal(total(pos) + total(neg), 1.0, 4)
	assert_array_almost_equal(total(values), 1.0, 4)

	positions, counts = entries(values, "chr1")
	assert_array_almost_equal(positions, [100, 200, 339, 419])
	assert_array_almost_equal(counts, [0.25, 0.125, 0.25, 0.125], 4)


def test_read_depth_with_scale_factor(stranded, bam, sizes):
	"""The documented way to make the bigWigs sum to a chosen value."""

	pos, neg = stranded(bam, "-s", sizes, "-r", "-sf", 1000000)

	assert_array_almost_equal(total(pos) + total(neg), 1000000, 4)

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 200])
	assert_array_almost_equal(counts, [250000, 125000], 4)


def test_read_depth_with_fragments_normalizes_to_one(stranded, bam, sizes):
	pos, neg = stranded(bam, "-s", sizes, "-f", "-r")

	assert_array_almost_equal(total(pos) + total(neg), 1.0, 4)


## Interval input: .bed, .bed.gz, .tsv, .tsv.gz


def test_bed_records_starts(stranded, bed, sizes):
	pos, _ = stranded(bed, "-s", sizes)

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 200])
	assert_array_almost_equal(counts, [2, 1], 4)

	positions, counts = entries(pos, "chr2")
	assert_array_almost_equal(positions, [50])
	assert_array_almost_equal(counts, [1], 4)


def test_bed_has_no_reverse_strand(stranded, bed, sizes):
	"""Interval files carry no strand, so everything lands on the plus file."""

	_, neg = stranded(bed, "-s", sizes)

	assert total(neg) == 0


def test_bed_fragments(stranded, bed, sizes):
	pos, _ = stranded(bed, "-s", sizes, "-f")

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 149, 200, 229])
	assert_array_almost_equal(counts, [2, 2, 1, 1], 4)

	positions, counts = entries(pos, "chr2")
	assert_array_almost_equal(positions, [50, 74])
	assert_array_almost_equal(counts, [1, 1], 4)


def test_bed_drops_chroms_not_in_sizes(stranded, bed, sizes):
	pos, _ = stranded(bed, "-s", sizes)

	assert "chrUN" not in pos
	assert total(pos) == 4


def test_bed_shifts(stranded, bed, sizes):
	pos, _ = stranded(bed, "-s", sizes, "-ps", 4, "-ns", -5, "-f")

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [104, 144, 204, 224])
	assert_array_almost_equal(counts, [2, 2, 1, 1], 4)


@pytest.mark.parametrize("other", ["bed_gz", "tsv", "tsv_gz"])
def test_interval_formats_agree(stranded, bed, sizes, other, request):
	"""The four interval extensions are the same parser and must agree."""

	other_path = request.getfixturevalue(other)

	from_bed, _ = stranded(bed, "-s", sizes, "-f", name="bed")
	from_other, _ = stranded(other_path, "-s", sizes, "-f", name="other")

	assert from_bed == from_other


def test_bed_float_coordinates(stranded, sizes, tmp_path):
	"""Fragment files exported from some pipelines carry float coordinates."""

	path = write_intervals(tmp_path / "floats.bed", [
		("chr1", "100.0", "150.0"),
		("chr1", "200.0", "230.0")
	])

	pos, _ = stranded(path, "-s", sizes, "-f")

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 149, 200, 229])
	assert_array_almost_equal(counts, [1, 1, 1, 1], 4)


def test_bed_extra_columns_are_ignored(stranded, sizes, tmp_path):
	path = write_intervals(tmp_path / "wide.bed", [
		("chr1", 100, 150, "peak1", 960, "+"),
		("chr1", 200, 230, "peak2", 12, "-")
	])

	pos, _ = stranded(path, "-s", sizes)

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 200])
	assert_array_almost_equal(counts, [1, 1], 4)


def test_empty_bed_produces_empty_bigwigs(run, sizes, tmp_path):
	path = write_intervals(tmp_path / "empty.bed", [])

	process = run(path, "-s", sizes)

	assert process.returncode == 0
	assert total(read_bigwig(tmp_path / "out.+.bw")) == 0


## Multiple input files and -p/--parallel


def test_two_files_are_concatenated(stranded, bam, sizes):
	pos, neg = stranded(bam, bam, "-s", sizes)

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 200])
	assert_array_almost_equal(counts, [4, 2], 4)

	positions, counts = entries(neg, "chr1")
	assert_array_almost_equal(positions, [339, 419])
	assert_array_almost_equal(counts, [4, 2], 4)


@pytest.mark.parametrize("n_jobs", [1, 2, 3])
def test_parallel_matches_serial(stranded, bam, sizes, n_jobs):
	pos, neg = stranded(bam, bam, bam, "-s", sizes, "-p", 1, name="serial")
	pos_p, neg_p = stranded(bam, bam, bam, "-s", sizes, "-p", n_jobs,
		name="parallel")

	assert pos == pos_p
	assert neg == neg_p


def test_parallel_unstranded_matches_serial(unstranded, bam, sizes):
	"""When unstranded, the two per-file dictionaries are the same object, so
	the merge across files has to avoid counting each file twice."""

	serial = unstranded(bam, bam, "-s", sizes, "-p", 1, name="serial")
	parallel = unstranded(bam, bam, "-s", sizes, "-p", 2, name="parallel")

	assert serial == parallel

	positions, counts = entries(serial, "chr1")
	assert_array_almost_equal(positions, [100, 200, 339, 419])
	assert_array_almost_equal(counts, [4, 2, 4, 2], 4)


def test_more_jobs_than_files(stranded, bam, sizes):
	pos, _ = stranded(bam, "-s", sizes, "-p", 4)

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 200])
	assert_array_almost_equal(counts, [2, 1], 4)


def test_file_order_does_not_matter(stranded, bam, bed, sizes):
	first, _ = stranded(bam, bed, "-s", sizes, name="first")
	second, _ = stranded(bed, bam, "-s", sizes, name="second")

	assert first == second


def test_mixed_bam_and_bed(stranded, bam, bed, sizes):
	pos, neg = stranded(bam, bed, "-s", sizes)

	# The BED entries land on the plus strand alongside the forward reads.
	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 200])
	assert_array_almost_equal(counts, [4, 2], 4)

	positions, counts = entries(neg, "chr1")
	assert_array_almost_equal(positions, [339, 419])
	assert_array_almost_equal(counts, [2, 1], 4)


def test_read_depth_spans_all_files(stranded, bam, sizes):
	pos, neg = stranded(bam, bam, "-s", sizes, "-r")

	assert_array_almost_equal(total(pos) + total(neg), 1.0, 4)


## -z/--zooms and -v/--verbose


@pytest.mark.parametrize("zooms", [0, 1, 4])
def test_zooms_do_not_change_values(stranded, bam, sizes, zooms):
	pos, neg = stranded(bam, "-s", sizes, name="default")
	pos_z, neg_z = stranded(bam, "-s", sizes, "-z", zooms, name="zoomed")

	assert pos == pos_z
	assert neg == neg_z


def test_verbose_reports_missing_chromosomes(run, bam, sizes):
	process = run(bam, "-s", sizes, "-v")

	assert process.returncode == 0
	assert "chrUN" in process.stdout + process.stderr


def test_quiet_by_default(run, bam, sizes):
	process = run(bam, "-s", sizes)

	assert process.returncode == 0
	assert "chrUN" not in process.stdout


## Combinations that occur in practice


def test_atac_bam_with_tn5_shift(stranded, bam, sizes):
	"""ATAC-seq reads, 5' cut sites, Tn5 offset applied per strand."""

	pos, neg = stranded(bam, "-s", sizes, "-ps", 4, "-ns", -5)

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [104, 204])
	assert_array_almost_equal(counts, [2, 1], 4)

	positions, counts = entries(neg, "chr1")
	assert_array_almost_equal(positions, [334, 414])
	assert_array_almost_equal(counts, [2, 1], 4)


def test_atac_fragments_unstranded(unstranded, tsv_gz, sizes):
	"""A 10x-style fragments.tsv.gz, both cut sites, one unstranded track."""

	values = unstranded(tsv_gz, "-s", sizes, "-f")

	positions, counts = entries(values, "chr1")
	assert_array_almost_equal(positions, [100, 149, 200, 229])
	assert_array_almost_equal(counts, [2, 2, 1, 1], 4)

	positions, counts = entries(values, "chr2")
	assert_array_almost_equal(positions, [50, 74])
	assert_array_almost_equal(counts, [1, 1], 4)


def test_atac_fragments_unstranded_shifted(unstranded, tsv_gz, sizes):
	values = unstranded(tsv_gz, "-s", sizes, "-f", "-ps", 4, "-ns", -5)

	positions, counts = entries(values, "chr1")
	assert_array_almost_equal(positions, [104, 144, 204, 224])
	assert_array_almost_equal(counts, [2, 2, 1, 1], 4)


def test_chip_seq_unstranded_cpm(unstranded, bam, sizes):
	"""ChIP-seq 5' ends, strand collapsed, normalized to counts per million."""

	values = unstranded(bam, "-s", sizes, "-r", "-sf", 1000000)

	assert_array_almost_equal(total(values), 1000000, 4)

	positions, counts = entries(values, "chr1")
	assert_array_almost_equal(positions, [100, 200, 339, 419])
	assert_array_almost_equal(counts, [250000, 125000, 250000, 125000], 4)


def test_merged_replicates_in_parallel(stranded, bam, sizes):
	"""Three replicate files pooled into one normalized pair of tracks."""

	pos, neg = stranded(bam, bam, bam, "-s", sizes, "-p", 3, "-r", "-sf",
		1000000)

	assert_array_almost_equal(total(pos) + total(neg), 1000000, 4)

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 200])
	assert_array_almost_equal(counts, [250000, 125000], 4)


def test_every_flag_together(unstranded, bam, sizes):
	"""Two files, both ends, shifted, pooled in parallel, CPM, zoomed."""

	values = unstranded(bam, bam, "-s", sizes, "-f", "-ps", 4, "-ns", -5,
		"-r", "-sf", 1000000, "-p", 2, "-z", 2)

	assert_array_almost_equal(total(values), 1000000, 4)

	positions, counts = entries(values, "chr1")
	assert_array_almost_equal(positions, [104, 144, 204, 224, 304, 334, 404,
		414])
	assert_array_almost_equal(counts, [125000, 125000, 62500, 62500, 125000,
		125000, 62500, 62500], 4)

	positions, counts = entries(values, "chr2")
	assert_array_almost_equal(positions, [54, 69, 84])
	assert_array_almost_equal(counts, [62500, 62500, 125000], 4)


@pytest.mark.parametrize("is_unstranded", [False, True])
@pytest.mark.parametrize("fragments", [False, True])
@pytest.mark.parametrize("shifts", [(0, 0), (4, -5)])
def test_total_signal_is_conserved(run, tmp_path, bam, sizes, is_unstranded,
	fragments, shifts):
	"""Across the whole flag cube, the amount of signal written depends only
	on how many reads were counted and whether both ends were recorded.

	Shifting moves positions and can collapse two ends onto one another, and
	going unstranded merges the two dictionaries, but neither may create or
	destroy signal.
	"""

	args = [bam, "-s", sizes, "-ps", shifts[0], "-ns", shifts[1]]
	if fragments:
		args.append("-f")
	if is_unstranded:
		args.append("-u")

	process = run(*args)
	assert process.returncode == 0, process.stderr

	if is_unstranded:
		written = total(read_bigwig(tmp_path / "out.bw"))
	else:
		written = (total(read_bigwig(tmp_path / "out.+.bw"))
			+ total(read_bigwig(tmp_path / "out.-.bw")))

	assert written == 8 * (2 if fragments else 1)


## Known bugs
#
# These describe how the tool should behave. They are skipped rather than
# asserting the current output, so that fixing the bug turns them green
# instead of requiring the test to be rewritten.


@pytest.mark.skip(reason="BUG: .sam passes the extension check but no branch "
	"in extract_reads handles it, so SAM input silently produces empty "
	"bigWigs instead of being read")
def test_sam_input_is_processed(stranded, sizes, tmp_path, bam):
	path = tmp_path / "test.sam"
	with pysam.AlignmentFile(str(bam), "rb") as infile:
		with pysam.AlignmentFile(str(path), "w", header=infile.header) as out:
			for alignment in infile:
				out.write(alignment)

	pos, neg = stranded(path, "-s", sizes)

	positions, counts = entries(pos, "chr1")
	assert_array_almost_equal(positions, [100, 200])
	assert_array_almost_equal(counts, [2, 1], 4)


@pytest.mark.skip(reason="BUG: a .sam path is never opened, so a nonexistent "
	"one exits 0 and writes empty bigWigs instead of erroring")
def test_missing_sam_file_errors(run, sizes, tmp_path):
	process = run(tmp_path / "does_not_exist.sam", "-s", sizes)

	assert process.returncode != 0
