#!/usr/bin/python
""" This program is designed to generate and validate primers from PAT-Seq for Re-PAT using primer3 and nesoni
Usage: get_primers.py gff_file	reference_file  gene_names_csv	output_file"""

import sys
import re
import glob
from nesoni import annotation, io, span_index
from subprocess import call 

def make_file_for_primer_3(gff_file, ref_file, names_file, output_file):

	gff_file = list(annotation.read_annotations(gff_file))
	print "\n Reading in the reference file"
	seq_dict = dict(io.read_sequences(ref_file))
	names_list = names_file
	config  = open("primer_config.txt").readlines()

	with open("regions_"+output_file, 'w') as out_f:
		for line in gff_file:
		    name = line.attr["Name"]
		    peak = line.attr["id"]
		    if name in names_list:
		        out_f.write ("SEQUENCE_ID="+ name + "_"+ peak + "\n")
		        out_f.write("SEQUENCE_TEMPLATE=" + line.shifted(-100, 0).get_seq(seq_dict) + "\n")
		        for cline in config:
		            out_f.write(cline.strip("\n")  + "\n")
		        out_f.write("="  + "\n")
       
def run_primer3(outfile):
	primer_3_input = "regions_"+ outfile
	print primer_3_input
	f = open ("primer_3_out_"+ outfile, "w")
	call (["primer3_core", primer_3_input ]) 
	#, stdout = f

def cat_reports (output_file):
	if len(glob.glob("./primer_3_reports")) == 0:
		call (["mkdir", "primer_3_reports"])
	if len(glob.glob("./*.for")) >0:
		call (["mv *.for primer_3_reports/"], shell =True)
	reports = glob.glob("./primer_3_reports/*")
	with open(output_file + "_potential_primer_list.csv", 'w') as out_f:
		out_f.write("Id" + "," + "Sequence" + "\n")
		for report in reports:
			print report 
			with open(report) as in_f:
				counter = 0
				for line in in_f:
					if counter < 4:
						counter += 1
						continue
					print line	
					split = line.strip("\n").split(" ")	
					if counter > 12:
						out_f.write (report.split("/")[-1] + "," + split [3] + "\n")
					else:
						out_f.write (report.split("/")[-1] + "," + split [4] + "\n")		      	
					counter += 1
					        
def csv_2_fa (output_file):
	with open(output_file + "test_primer_list.fa", "w") as out_f:
	    with open(output_file + "_potential_primer_list.csv") as in_f:
	      counter = 0
	      for line in in_f:
	        if counter > 0:
	          split = line.strip("\n").split(",")
	          print split
	          name = split[0]
	          seq = split[1]
	          out_f.write(">"+name + "\n"+ seq + "\n")

	        counter = counter +1


def interrogate_exonerate(output_file):
	with open ("temp_exon_rep") as in_f:
		hit_count = 0
		for line in in_f:
			print line
			if "C4 Alignment:" in line:
				hit_count +=1
				print hit_count
	if hit_count == 0:
		print "Warning! Primer not detected in reference sequence"
		return False
	elif hit_count == 1:
		return True
	else:
		return False


def xoneratebby (output_file, reference_file):
	with open(output_file + "final_primer_list.csv","w") as final:
		with open(output_file + "test_primer_list.fa") as in_f:
			for line in in_f:			
				if line[0] == ">":
					name = line.strip(">")
					if "replace" in locals() and name == replace:
						continue
					seq = in_f.next()
					print name, seq
					with open("temp.fa", "w") as ex_f:
						ex_f.write(">" + name + seq)
					f = open ("temp_exon_rep", "w")
					call (["exonerate", "temp.fa", reference_file], stdout = f)
					if interrogate_exonerate(output_file) == True:
						final.write(name.strip("\n") + "," + seq.strip("\n") + "\n" )
						replace = name
					else:
						continue

					#"--model", "affine:local", "--exhaustive", "yes", "--score", "115",
				#elif line != name:
					#name = line
					#seq = line.next

				








					

gff = sys.argv[1]
ref = sys.argv[2]
names = sys.argv[3]
output_file = sys.argv[4]

make_file_for_primer_3(gff, ref, names, output_file)
run_primer3(output_file)
cat_reports(output_file)
csv_2_fa(output_file)
#xoneratebby (output_file, ref)
# you might need to deal with the rev strand. 

def generate_report(outfile):
	with open ("primer_3_out_" + outfile) as in_f:
		with open ("report_" + outfile, "w") as out_f:
			out_f.write("Id" + "\t" + "Sequence"+ "\t" + "TM" + "\t" + "GC percent" + "\t" + "Hairpin TH" + "\t" + "End TH" + "\n")
			counter = 0
			for line in in_f:							
				split = line.strip("\n").split("=")
				if counter % 6 == 0:
					out_f.write("\n")
				elif re.match(split[0], "SEQUENCE_ID"):
					out_f.write(line.split("=")[1] + "\t")
				elif re.match(split[0], "*_SEQUENCE"):
					print split[1]
					out_f.write(split[1] + "\t")
				elif re.match(split[0], "*_TM"):
					out_f.write(split[1] + "\t")
				elif re.match(split[0], "*_GC_PERCENT"):
					out_f.write(split[1] + "\t") 
				elif re.match(split[0], "*_HAIRPIN_TH"):
					print ('try')
					out_f.write(split[1] + "\t")
				elif re.match(split[0], "*_SELF_END_TH"):
					print "tu"
					out_f.write(split[1] + "\t") 
				
				counter += 1