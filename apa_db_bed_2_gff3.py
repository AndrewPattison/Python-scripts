#!/usr/bin/env python

"""Usage python apa_db_Bed_2_gff.py <in.bed> <out.gff>"""
import sys
def reorder_dp (input_file, outfile):
    prog_name = "bed_reorderer"
    dot = "."
    with open(outfile, 'w') as out_f:
            out_f.write("##gff-version 3 \n")
            with open(input_file) as f:
                for line in f:
                    split_line = line.strip().split("\t")
                    chromo = split_line[0]
                    strand = split_line[5]
                    start = str(int (split_line[1])+1)
                    end = split_line[2]
                    score = "."
                    name = split_line[3]
                    mirs= split_line[-1]
                    supporting_reads = split_line[4]
                    if strand == "-":
                        colour = "#0000ff"

                    else:
                        colour = "#00ff00"   

                    out_line = chromo +"\t" +prog_name +"\t"+ "apa_db_Bed_2_gff" + "\t"+ start+ "\t"+ end+ "\t"+ score  +"\t"+strand+"\t"+ "."+"\t"+ "color=" +colour+";"+ "name="+ name +";"+"APADB_supporting_reads="+ supporting_reads+";"+"lost_mirs ="+mirs+"\n"
                    out_f.write(out_line)                       
                               


if __name__ == "__main__":

    in_file = sys.argv[1]

    out_fie= sys.argv[2]

    reorder_dp (in_file,out_fie )