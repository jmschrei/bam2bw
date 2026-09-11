# synthetic.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com>

import gzip

import pysam


# Every input that `bam2bw` accepts is built here from an explicit, tiny read
# layout rather than being committed as a binary fixture. The point is that the
# expected bigWig values in the test file can be derived by hand from the
# layout: a read that starts at 100 and is 50bp long contributes its 5' end at
# 100 and its 3' end at 149, and nothing in between. Committing a real BAM
# would make those numbers unauditable.


def write_chrom_sizes(path, chrom_sizes):
	"""Write a two-column chrom_sizes file.

	Parameters
	----------
	path: str or pathlib.Path
		The file to write to.

	chrom_sizes: list of (str, int) tuples
		The chromosome names and their lengths, in the order they should
		appear. `bam2bw` preserves this order in the bigWig header.

	Returns
	-------
	path: str or pathlib.Path
		The path that was written to, for convenience.
	"""

	with open(path, "w") as outfile:
		for chrom, size in chrom_sizes:
			outfile.write("{}\t{}\n".format(chrom, size))

	return path


def write_fasta(path, chrom_sizes, line_width=60):
	"""Write a FASTA file whose sequences have the requested lengths.

	`bam2bw` accepts a FASTA in place of a chrom_sizes file and reads the
	lengths out of it with pyfaidx, so only the lengths matter here. The
	sequence content is a repeating ACGT so that the file looks like real
	DNA to anything that inspects it.

	Note that pyfaidx writes a `.fai` index next to this file the first time
	it is opened, so the path should be inside a scratch directory.

	Parameters
	----------
	path: str or pathlib.Path
		The file to write to.

	chrom_sizes: list of (str, int) tuples
		The chromosome names and the length of the sequence to write for each.

	line_width: int, optional
		The number of bases per line of sequence. Default is 60.

	Returns
	-------
	path: str or pathlib.Path
		The path that was written to, for convenience.
	"""

	with open(path, "w") as outfile:
		for chrom, size in chrom_sizes:
			sequence = ("ACGT" * (size // 4 + 1))[:size]

			outfile.write(">{}\n".format(chrom))
			for i in range(0, size, line_width):
				outfile.write(sequence[i:i+line_width] + "\n")

	return path


def write_bam(path, chrom_sizes, reads, unmapped=()):
	"""Write a BAM file containing the given reads.

	Each read is placed with an exact-match CIGAR so that its reference span
	is precisely `start` to `start + length`, which is what makes the expected
	5' and 3' positions hand-computable. The header can list more chromosomes
	than the chrom_sizes file passed to `bam2bw` does, which is how the
	missing-chromosome filter gets exercised.

	Parameters
	----------
	path: str or pathlib.Path
		The file to write to.

	chrom_sizes: list of (str, int) tuples
		The chromosomes to put in the BAM header. Reads may only refer to
		chromosomes in this list.

	reads: list of (str, int, int, bool) tuples
		One tuple per read, giving the chromosome, the 0-based reference
		start, the length of the alignment, and whether the read is on the
		reverse strand.

	unmapped: list of (str, int, int) tuples, optional
		One tuple per unmapped read to append, giving a chromosome, a start,
		and a length. The reads carry the unmapped flag but are still placed,
		the way a read whose mate aligned is placed next to its mate in a real
		BAM. Placing them matters for the tests: an unmapped read with no
		reference would be dropped by the missing-chromosome filter before the
		unmapped filter ever saw it, which would make the two filters
		impossible to tell apart. Default is an empty list.

	Returns
	-------
	path: str or pathlib.Path
		The path that was written to, for convenience.
	"""

	header = {
		'HD': {'VN': '1.6', 'SO': 'unsorted'},
		'SQ': [{'SN': chrom, 'LN': size} for chrom, size in chrom_sizes]
	}

	tids = {chrom: i for i, (chrom, _) in enumerate(chrom_sizes)}

	with pysam.AlignmentFile(str(path), "wb", header=header) as outfile:
		for i, (chrom, start, length, is_reverse) in enumerate(reads):
			alignment = pysam.AlignedSegment(outfile.header)
			alignment.query_name = "read{}".format(i)
			alignment.query_sequence = "A" * length
			alignment.query_qualities = pysam.qualitystring_to_array("I" * length)
			alignment.flag = 0
			alignment.reference_id = tids[chrom]
			alignment.reference_start = start
			alignment.mapping_quality = 60
			alignment.cigar = ((0, length),)
			alignment.is_reverse = is_reverse

			outfile.write(alignment)

		for i, (chrom, start, length) in enumerate(unmapped):
			alignment = pysam.AlignedSegment(outfile.header)
			alignment.query_name = "unmapped{}".format(i)
			alignment.query_sequence = "A" * length
			alignment.query_qualities = pysam.qualitystring_to_array("I" * length)
			alignment.flag = 4
			alignment.reference_id = tids[chrom]
			alignment.reference_start = start
			alignment.cigar = ((0, length),)

			outfile.write(alignment)

	return path


def write_intervals(path, entries):
	"""Write a BED/tsv file of intervals, gzipping it if the path says to.

	The fields are written through `str` exactly as given, so a test can pass
	floats or strings to exercise the coordinate parsing, and can pass tuples
	longer than three to exercise the columns that `bam2bw` ignores.

	Parameters
	----------
	path: str or pathlib.Path
		The file to write to. If it ends in `.gz` the file is gzipped.

	entries: list of tuples
		One tuple per line. The first three fields are the chromosome, the
		start, and the end; any further fields are written but ignored by
		`bam2bw`.

	Returns
	-------
	path: str or pathlib.Path
		The path that was written to, for convenience.
	"""

	if str(path).endswith(".gz"):
		outfile = gzip.open(path, "wt")
	else:
		outfile = open(path, "w")

	with outfile:
		for entry in entries:
			outfile.write("\t".join(str(field) for field in entry) + "\n")

	return path
