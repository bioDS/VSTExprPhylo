#!/usr/bin/env python3

import sys

xml_template_filename=sys.argv[1]
taxonset_filename=sys.argv[2]
keys_filename=sys.argv[3]
counts_filename=sys.argv[4]
basename=sys.argv[5]
minor_dimension=sys.argv[6]
output_filename=sys.argv[7]

with open(keys_filename, 'r') as keys_file:
    keys=keys_file.readline().strip()
    
with open(xml_template_filename, 'r') as f:
    xml_lines = f.readlines()

with open(output_filename, 'w') as out_file:
    for line in xml_lines:
        if "insert_taxonset" in line:
            with open(taxonset_filename, 'r') as taxonset_file:
                for taxon in taxonset_file.readlines():
                    out_file.write(taxon)
        elif "insert_keys" in line:
            new_line=line.replace("insert_keys", keys).replace("insert_minor_dimension", minor_dimension)
            out_file.write(new_line)
        elif "insert_counts" in line:
            with open(counts_filename, 'r') as counts_file:
                counts_line=counts_file.readline()
            out_file.write(counts_line)
        elif "filename" in line:
            new_line=line.replace("filename", basename)
            out_file.write(new_line)
        else:
            out_file.write(line)
        
        


