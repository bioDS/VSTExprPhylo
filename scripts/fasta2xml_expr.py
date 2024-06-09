#!/usr/bin/env python3

import sys
from Bio import SeqIO

xml_template_file=sys.argv[1]
expr_fasta_file=sys.argv[2]
output_file=sys.argv[3]

expr_fasta_sequences = SeqIO.parse(open(expr_fasta_file),'fasta')

with open(xml_template_file, 'r') as f:
    xml_lines = f.readlines()

with open(output_file, 'w') as out_file:
    for line in xml_lines:
        if "{{expr_seq}}" in line:
            for fasta in expr_fasta_sequences:
                name, sequence = fasta.id, str(fasta.seq)
                out_file.write("<sequence taxon=\"" + name + "\">\n")
                out_file.write(sequence+"\n")
                out_file.write("</sequence>\n")
        else:
            out_file.write(line)
        
        


