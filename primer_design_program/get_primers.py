#!/usr/bin/python
""" This program is designed to generate and validate primers from PAT-Seq for Re-PAT using primer3 and nesoni
Usage: get_primers.py gff_file  reference_file  gene_names_text_file (line separated) output_file"""

import sys
import re
import glob
from nesoni import annotation, io, span_index
from subprocess import call 


# This function uses Paul Harrison's annotation class to find regions to give them to primer3
# Regions are from input genes
# They are written alongside parameters set in the primer3 config file 
def make_file_for_primer_3 (gff_file, ref_file, names_file, output_file):

    gff_file = list(annotation.read_annotations(gff_file))
    print "\n Reading in the reference file"
    seq_dict = dict(io.read_sequences(ref_file))
    names_file = open(names_file).readlines()
    config  = open("primer_config.txt").readlines()

    with open("regions_" + output_file, 'w') as out_f:
        for name in names_file:
            sname = name.strip("\n")
            found = False
            for line in gff_file:
                gff_name = line.attr.get ("Name", "No_name")
                peak = line.attr.get ("id", "No_id")
                if sname in gff_name.split("/"):
                    out_f.write ("SEQUENCE_ID="+ gff_name + "_" + peak + "\n")
                    out_f.write("SEQUENCE_TEMPLATE=" + line.shifted(-100, 0).get_seq(seq_dict) + "\n")
                    found = True
                    for cline in config:
                        out_f.write(cline.strip("\n")  + "\n")
                    out_f.write("="  + "\n")
            if found ==False:
                print "Could not find the gene " + sname + " in the gff file"

# Runs primer3 on the list of sequences and leaves and output for each one
# in the current working directory (primer3 does this)
def run_primer3(outfile):
    primer_3_input = "regions_"+ outfile
    f = open ("primer_3_out_"+ outfile, "w")
    call (["primer3_core", primer_3_input]) 

# Grabs all the reports in the current working directory and throws them into a new one
# which it makes if it doesn't exist
def cat_reports (output_file):
    if len(glob.glob("./primer_3_reports")) == 0:
        call (["mkdir", "primer_3_reports"])
    if len(glob.glob("./*.for")) >0:
        call (["rm ./primer_3_reports/*.for"], shell =True)
        call (["mv *.for primer_3_reports/"], shell =True)
    reports = glob.glob("./primer_3_reports/*")
    with open(output_file + "_potential_primer_list.csv", 'w') as out_f:
        out_f.write("Id" + "," + "Sequence" + "\n")
        for report in reports:
            with open(report) as in_f:
                nlines = in_f.readlines()
                print "Primer3 has found " + str(len(nlines)) + " potential primers in " + report.strip(".for").split("/")[-1]
                counter = 0
                for line in nlines:
                    if counter < 4:
                        counter += 1
                        continue
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
              name = split[0]
              seq = split[1]
              if len(seq) >1:
                out_f.write(">"+name + "\n"+ seq + "\n")

            counter = counter +1


def interrogate_exonerate(output_file):
    hit_count = 0
    with open ("temp_exon_rep") as in_f:
        for line in in_f:
            print line
            if "C4 Alignment:" in line:
                hit_count +=1
            if hit_count >1:
                print "More than one region found"
                return False
    if hit_count == 1:
        print "Unique primer found"
        return True
    else:
        print "Warning! Primer not detected in reference sequence"
        return False
        


def xoneratebby (output_file, reference_file):
    with open(output_file + "final_primer_list.csv","w") as final:
        final.write("Id" + "," + "Primer" + "\n")
        with open(output_file + "test_primer_list.fa") as in_f:
            counter = 0 
            found = False 
            for line in in_f:         
                if line[0] == ">":
                    name = line.strip(">")
                    if found == True and name == last_name:
                        continue
                    if counter > 0 and found == False and name != last_name:
                        print "Warning! No suitable regions found for " + name  + "\n"
                        final.write(name.strip(".for\n") + "," + "No_seq_found" + "\n" )
                    seq = in_f.next()
                    print name, seq
                    with open("temp.fa", "w") as ex_f:
                        ex_f.write(">" + name + seq)
                    f = open ("temp_exon_rep", "w")
                    call (["exonerate","-n", "10", "--percent", "95", "-m", "affine:local" , "temp.fa", reference_file], stdout = f)
                    if interrogate_exonerate(output_file) == True:
                        if counter > 0 and seq == last_seq:
                            print "Same sequence found for two close primers"
                            final.write(name.strip(".for\n") + "," + "Same_as_previous" + "\n")                            
                        else:
                            final.write(name.strip(".for\n") + "," + seq.strip("\n") + "\n" )
                            
                        found = True
                        last_name = name
                        last_seq =seq 

                    else:
                        found = False
                        last_name = name
                        counter = 0
                        last_seq =seq 
                        continue
                    counter += 1

                    #"--model", "affine:local", "--exhaustive", "yes", "--score", "115",
                    #tail-tools primer-gff prefix  /data/reference/tail-tools/h_sapiens_ensembl_82/ tes_5final_primer_list.csv  

def primer_gff(output_file, ref):
    call (["tail-tools", "primer-gff", output_file,  ref.strip("reference.fa"),  output_file + "final_primer_list.csv"])

gff = sys.argv[1]
ref = sys.argv[2]
names = sys.argv[3]
output_file = sys.argv[4]

make_file_for_primer_3(gff, ref, names, output_file)
run_primer3(output_file)
cat_reports(output_file)
csv_2_fa(output_file)
xoneratebby (output_file, ref)
primer_gff(output_file, ref)
#you might need to deal with the rev strand. Might not
# Build in a warning if primer3 finds nothing for a site. 
# Better warnings like a sys.stderror 
# More checks of the data
# Maybe more options
# Show how may primers primer3 cam up with.
# Turn off blast low complexity regions 