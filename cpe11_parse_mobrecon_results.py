#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <user@localhost>
Date   : 2024-10-24
Purpose: Parse mobrecon data.
"""

import argparse


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Parse mobrecon data.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-c',
                        '--mobrecon_contigs',
                        help='mobrecon contigs results file',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-m',
                        '--mobtyper_results',
                        help='mobrecon mobtyper results file',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-s',
                        '--sample',
                        help='sample/isolate name',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-i',
                        '--input_file',
                        help='output file from cpe11_get_carba_contigs.py',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-a',
                        '--assembler',
                        help='unicycler or skesa, since they name contigs differently',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-oc',
                        '--output_contigs',
                        help='carbapenemase containing contigs output',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-op',
                        '--output_plasmids',
                        help='carbapenemase containing plasmids output',
                        metavar='str',
                        type=str,
                        required=True)

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()

    plasmids = []
    contigs = []

    # Parse carbapenemase contigs:
    carba_contigs = []

    with open(args.input_file, "r") as infile1:
        for line in infile1:
            if args.sample in line:
                line_elements = line.strip().split("\t")
                carba_contigs.append(line_elements[1])
    
    mobrecon_contigs = []
    with open(args.mobrecon_contigs, "r") as infile2:
        for line in infile2:
            mobrecon_contigs.append(line.strip())

    if args.assembler == "unicycler":
        for contig in carba_contigs:
            for line in mobrecon_contigs:
                if contig.split("_")[1] in line.split("\t")[4].split(" ")[0]:
                    contig_line = line.strip()
                    contigs.append(contig_line)
                    plasmid = line.split("\t")[2]

    elif args.assembler == "skesa":
        for contig in carba_contigs:
            for line in mobrecon_contigs:
                # Contig_103_57.6422 format for skesa from mobrecon (for some reason)
                if contig.split("_")[1] in line.split("\t")[4].split("_")[1]:
                    contig_line = line.strip()
                    contigs.append(contig_line)
                    plasmid = line.split("\t")[2]

    if plasmid != "-":
        with open(args.mobtyper_results, "r") as infile3:
            for line in infile3:
                if plasmid in line:
                    plasmid_line = line.strip()
                    plasmids.append(plasmid_line)
    
    with open(args.output_contigs, "a") as outfile1:
        to_write = "\n"+"\n".join(contigs)
        outfile1.write(to_write)
    
    with open(args.output_plasmids, "a") as outfile2:
        to_write2 = "\n"+"\n".join(plasmids)
        outfile2.write(to_write2)






# --------------------------------------------------
if __name__ == '__main__':
    main()
