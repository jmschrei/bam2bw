from setuptools import setup

setup(
    name='bam2bw',
    version='0.3.1',
    author='Jacob Schreiber',
    author_email='jmschreiber91@gmail.com',
    scripts=['bam2bw'],
    url='https://github.com/jmschrei/bam2bw',
    license='LICENSE.txt',
    description='A command-line tool for converting SAM/BAM files into stranded bp resolution BigWig files. ',
    install_requires=[
        "numpy",
        "pysam",
		"pyfaidx",
        "pyBigWig",
        "tqdm"
    ],
)