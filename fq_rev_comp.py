#!/usr/bin/env python
from __future__ import print_function
from Bio.Seq import Seq
import sys


"""Tool for getting the reverse complement of a fastq file. By Andrew Pattison
Usage: python fq_rev_comp.py Input_filename Output_filename"""

# If I use this again iterate over the file one line at a time to save memory. 
def rev_comp (fqfile, output_filename):
	with open(output_filename, 'w') as f:

		for line in xrange(0,len(fqfile)):
			if line % 4 == 1:
				assert fqfile[line][0] in "ATCGN", "Either we are on the wrong line or there is a base other than ATCG or N"	
				no_re_line = fqfile[line].rstrip()
				new_line = Seq(no_re_line)
				rev_line = new_line.reverse_complement()
				f.write(str(rev_line)+"\n")

			elif line % 4 == 3:
				new_qual_line = fqfile[line].rstrip()
				new_qual_line = new_qual_line[::-1]
				f.write(new_qual_line+"\n")
			else:
				f.write(fqfile[line])


if __name__ == "__main__":

	fastq_file = sys.argv[1]

	output_filename = sys.argv[2]

	fqfile = open (fastq_file, "rU").readlines()

	rev_comp(fqfile, output_filename)