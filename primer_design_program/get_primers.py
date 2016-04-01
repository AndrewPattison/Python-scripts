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
    # check for a tmp direcory 
    if len(glob.glob("./tmp")) == 0:
        call (["mkdir", "tmp"])

    gff_file = list(annotation.read_annotations(gff_file))
    print "\n Reading in the reference file"
    seq_dict = dict(io.read_sequences(ref_file))
    names_file = open(names_file).readlines()
    config  = open("primer_config.txt").readlines()

    with open("tmp/regions_" + output_file, 'w') as out_f:
        for name in names_file:
            sname = name.strip("\n")
            found = False
            for line in gff_file:
                gff_name = line.attr.get ("Name", "No_name")
                peak = line.attr.get ("id", "No_id")
                if sname in gff_name.split("/"):
                    out_f.write ("SEQUENCE_ID="+ gff_name.replace("/", "_") + "_" + peak + "\n")
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
    primer_3_input = "tmp/regions_"+ outfile
    f = open ("tmp/primer_3_out_"+ outfile, "w")
    call (["primer3_core", primer_3_input], stdout = f) 

# Grabs all the reports in the current working directory and throws them into a new one
# which it makes if it doesn't exist
def cat_reports (output_file):
    if len(glob.glob("./primer_3_reports")) == 0:
        call (["mkdir", "primer_3_reports"])
    if len(glob.glob("./*.for")) >0:
        call (["rm ./primer_3_reports/*.for"], shell =True)
        call (["mv *.for primer_3_reports/"], shell =True)
    reports = glob.glob("./primer_3_reports/*")
    with open("tmp/" + output_file + "_potential_primer_list.csv", 'w') as out_f:
        out_f.write("Id" + "," + "Sequence" + "\n")
        for report in reports:
            with open(report) as in_f:
                nlines = in_f.readlines()
                with open ("tmp/report.txt", "a") as report_txt:
                    report_txt.write("Primer3 has found " + str(len(nlines)) + " potential primers in " + report.strip(".for").split("/")[-1])
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

def check_sequence_comp(seq):
    with open ("tmp/report.txt", "w") as report_txt:
        total = len(seq)
        a = float(seq.count("A"))/total
        t = float(seq.count("T"))/total
        g = float(seq.count("G"))/total
        c = float(seq.count("C"))/total
        if a or t or c or g < 0.15:
            report.write("Sequence " + seq + " removed due to poor sequence composition" + "\n")
            return False
        else:
            return True

def check_short_seq_comp(seq):
    with open ("tmp/report.txt", "a") as report_txt:
        short_seq = seq[:-16]
        a = seq.count("A")
        t = seq.count("T")
        g = seq.count("G")
        c = seq.count("C")

        if a or t or c or g < 2:
            report.write("Sequence " + seq + " removed due to poor end of primer sequence composition" + "\n")
            return False
        else:
            return True

                            
def csv_2_fa (output_file):
    with open("tmp/"+ output_file + "test_primer_list.fa", "w") as out_f:
        with open("tmp/" + output_file + "_potential_primer_list.csv") as in_f:
          counter = 0
          removed_counter = 0 
          for line in in_f:
            if counter > 0:
              split = line.strip("\n").split(",")
              name = split[0]
              seq = split[1]
              if len(seq) >1:
                if check_sequence_comp(seq) and check_short_seq_comp(seq):
                    out_f.write(">"+name + "\n"+ seq + "\n")
                removed_counter += 1

            counter = counter +1
        print removed_counter + " total primers removed due to unequal base composition"


def interrogate_exonerate():
    with open ("tmp/report.txt", "a") as report:
        hit_count = 0
        with open ("tmp/temp_exon_rep") as in_f:
            for line in in_f:
                if "Sbjct" in line:
                    hit_count +=1
                if hit_count >1:
                    report.write( "More than one region found")
                    return False
        
        if hit_count == 1:
            report.write("Unique primer found")
            return True
        else:
            report.write("Warning! Primer not detected in reference sequence")
            return False

def interrogate_short():
    hit_count = 0
    with open("tmp/short_exon_rep") as in_f:
        for line in in_f:
            if "Sbjct" in line:
                hit_count +=1
    print str(hit_count) + "hits found"
    return(hit_count)

def run_blast(name, seq, reference_file):
    with open("tmp/temp.fa", "w") as ex_f:
        ex_f.write(">" + name + "\n" + seq)
    f = open ("tmp/temp_exon_rep", "w")
    call(["blastn", "-task", "blastn-short", "-word_size", str(len(seq)-1), "-num_threads", "4", "-dust","no", "-perc_identity", "100", "-query", "tmp/temp.fa", "-db", reference_file], stdout= f)

def check_short_seq(name,seq, ref):
    print len(seq [-13:])
    with open("tmp/temp_short.fa", "w") as short_f:
        short_f.write(">" + name  + "\n" + seq.strip("\n") [-15:])
    f = open ("tmp/short_exon_rep", "w")
    call(["blastn", "-task", "blastn-short", "-word_size", "11", "-num_threads", "4", "-dust","no", "-perc_identity", "100", "-query", "tmp/temp_short.fa", "-db", ref], stdout= f)

        

def xoneratebby (output_file, reference_file):
    with open("tmp/" + output_file + "final_primer_list.csv","w") as final:
        final.write("Id" + "," + "Primer" + "\n")
        with open("tmp/" + output_file + "test_primer_list.fa") as in_f:
            for line in in_f:         
                if line[0] == ">":
                    name = line.strip(">")
                    seq = in_f.next()
                    run_blast(name, seq, reference_file)
                    if interrogate_exonerate() == True:
                        final.write (name.strip(".for\n") + "," + seq)

def best_primer_by_primer3(output_file):
    with open(output_file + "_primer3_suggested.csv", "w") as out_f:
        out_f.write("Id,Primer" + "\n")
        with open("tmp/" + output_file + "final_primer_list.csv") as final:
            count = 0
            last_name = ""
            for line in final:
                if line[0] != "\n":
                    split = line.split(",")
                    name = split[0]
                    if count > 0:
                        if name in last_name:
                            continue
                        else:
                            out_f.write(line)
                    last_name = name
                    count += 1

def best_by_short_seq(output_file, ref):
    with open(output_file + "_blast_suggested.csv", "w") as out_f:
        out_f.write("Id,Primer" + "\n")
        with open("tmp/" + output_file + "final_primer_list.csv") as final:
            last_name = ""
            best = ""
            name =""
            count = 0
            min_hits = 1000000
            for line in final:
                if count == 0:
                    count += 1
                    continue
                if line[0] != "\n":
                    split = line.split(",")
                    if name != last_name:
                        out_f.write(best)
                    name = split[0]
                    seq = split[1].strip("\n")
                    print name ,seq
                    check_short_seq(name,seq, ref)
                    hits = interrogate_short()
                    if hits < min_hits:
                        min_hits = hits
                        best = line

                    last_name = name
                    count += 1

def primer_gff(output_file, ref):
    call (["tail-tools", "primer-gff", output_file + "_primer3_suggested",  ref.strip("reference.fa"),  output_file + "_blast_suggested.csv"])

def primer_gff_2(output_file, ref):
    call (["tail-tools", "primer-gff", output_file + "_BLAST_suggested",  ref.strip("reference.fa"),  output_file + "_.csv"])

gff = sys.argv[1]
ref = sys.argv[2]
names = sys.argv[3]
output_file = sys.argv[4]

make_file_for_primer_3(gff, ref, names, output_file)
run_primer3(output_file)
cat_reports(output_file)
csv_2_fa(output_file)
xoneratebby (output_file, ref)
best_primer_by_primer3(output_file)
best_by_short_seq(output_file, ref)
primer_gff(output_file, ref)
primer_gff2(output_file, ref)
# Build in a warning if primer3 finds nothing for a site. 
# Better warnings like a sys.stderror 
# More checks of the data
# Maybe more options
# Show how may primers primer3 cam up with.
