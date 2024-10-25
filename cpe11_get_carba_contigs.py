#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <user@localhost>
Date   : 2024-10-24
Purpose: Parse mobrecon results.
"""

import argparse


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Parse mobrecon results.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    
    parser.add_argument('-c',
                        '--card_input',
                        help='card input text file',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-g',
                        '--gff',
                        help='isolate gff3 file from bakta',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-s',
                        '--sample',
                        help='sample/isolate name',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-o',
                        '--outfile',
                        help='output file; will be appended to',
                        metavar='str',
                        type=str,
                        required=True)

    return parser.parse_args()


# --------------------------------------------------


def parse_rgi(card_handle, bakta_gff_handle, isolate, output_file):

    gene_identifiers = []
    genes = []
    output_lines = []
    
    with open(card_handle, "r") as infile1:
        for line in infile1:
            line_elements = line.strip().split("\t")
            if any(a in line_elements[8] for a in ["IMI","IMP","KPC","VIM","NDM","NMC","OXA","GES"]):
                gene_identifier = line_elements[1].split("_")[0]+"_"+line_elements[1].split("_")[1]
                gene = line_elements[8]

                gene_identifiers.append(gene_identifier)
                genes.append(gene)
    
    # if gene_identifier != "":
    for gene_identifier, gene in zip(gene_identifiers, genes):
    
        with open(bakta_gff_handle, "r") as infile2:
            for line in infile2:
                if gene_identifier in line:
                    line_elements = line.strip().split("\t")
                    if line_elements[2] == "gene":
                        contig = line_elements[0]
                        output_line = f"{isolate}\t{contig}\t{gene_identifier}\t{gene}"
                        output_lines.append(output_line)

    to_write = "\n" + "\n".join(output_lines)
    with open(output_file, "a") as outfile1:
        outfile1.write(to_write)


def main():
    """Generate file indicating contig of carbapenemase (if present in assembly)"""

    args = get_args()

    parse_rgi(args.card_input, args.gff, args.sample, args.outfile)



# --------------------------------------------------
if __name__ == '__main__':
    main()
