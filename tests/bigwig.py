# bigwig.py
# Contact: Jacob Schreiber <jmschreiber91@gmail.com>

import pyBigWig


# `bam2bw` writes one entry per non-zero position with span=1, so reading a
# bigWig back gives a sparse mapping of position to value. Everything the tests
# assert on is expressed through that mapping rather than through the file
# bytes, because the byte layout is pyBigWig's business and can change between
# versions without the counts being wrong.


def read_bigwig(path):
	"""Read a bigWig into a nested dictionary of positions and values.

	Every chromosome in the header appears as a key even when it carries no
	entries, so that a test can assert a chromosome is empty rather than
	having to assert it is absent.

	Parameters
	----------
	path: str or pathlib.Path
		The bigWig file to read.

	Returns
	-------
	values: dict of str to dict of int to float
		A mapping from chromosome name, to a mapping from position to the
		value recorded there. Positions with no signal are not present.
	"""

	values = {}

	bw = pyBigWig.open(str(path))
	for chrom in bw.chroms():
		values[chrom] = {}

		for start, end, value in (bw.intervals(chrom) or ()):
			for position in range(start, end):
				values[chrom][position] = value

	bw.close()

	return values


def entries(values, chrom):
	"""Split one chromosome of a read bigWig into sorted positions and values.

	This exists so that assertions can be written as two flat lists compared
	with `assert_array_almost_equal`, which is both easier to read than a dict
	literal and reports the first differing position when it fails.

	Parameters
	----------
	values: dict of str to dict of int to float
		The return value of `read_bigwig`.

	chrom: str
		The chromosome to pull out.

	Returns
	-------
	positions: list of int
		The positions carrying signal, in ascending order.

	counts: list of float
		The value at each of those positions, in the same order.
	"""

	positions = sorted(values[chrom])
	counts = [values[chrom][position] for position in positions]

	return positions, counts


def total(values):
	"""Sum every value across every chromosome of a read bigWig.

	Parameters
	----------
	values: dict of str to dict of int to float
		The return value of `read_bigwig`.

	Returns
	-------
	total: float
		The sum of all recorded values.
	"""

	return sum(sum(chrom.values()) for chrom in values.values())


def chrom_lengths(path):
	"""Read the chromosome lengths out of a bigWig header.

	The lengths never appear in the values, so a sizes file or FASTA that
	gives the wrong length produces a correct-looking track with a wrong
	header. This is the only way to see that.

	Parameters
	----------
	path: str or pathlib.Path
		The bigWig file to read.

	Returns
	-------
	lengths: dict of str to int
		A mapping from chromosome name to the length in the header.
	"""

	bw = pyBigWig.open(str(path))
	lengths = dict(bw.chroms())
	bw.close()

	return lengths
