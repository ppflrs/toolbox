#!/usr/bin/env python

import sys
import argparse

arg_parser = argparse.ArgumentParser(description='Linearize a FASTA file and print to the standard output.')
arg_parser.add_argument("--tab", "-t", action="store_true", help='Print sequences in tab separated format.')
arg_parser.add_argument("--len", "-l", action="store_true", help='Print lenght of sequence after the sequence header. Only for tab separated format.')
arg_parser.add_argument("--rename", "-r", type=str, help='<tsv_file> Rename the headers to "seq_1, seq_2, etc." format. The original headers are stored in tsv_file.')
arg_parser.add_argument("input_file",type=str, help='<input_file> full path to fasta file.')

args = arg_parser.parse_args()

infile = args.input_file
tab_print = args.tab
len_print = args.len
rename_seqs_file = args.rename

if rename_seqs_file != None:
    rename_seqs = True
    seq_counter = 0
    if rename_seqs_file == infile:
        sys.exit("Warning: This will delete the input file, please change the renamed headers mapping file.")
else:
    rename_seqs = False

if len_print == True:
    if tab_print == False:
        sys.exit("Sorry --len option only valid if used with --tab.")

with open(infile, 'rt') as fh:
    entry = ''
    seq = ''
    header = ''
    if rename_seqs == True:
        new_header = ''
        with open(rename_seqs_file, "wt") as renamed_seqs_fh:
            if tab_print == False:
                    for line in fh:
                            if line[0] == '>':
                                    original_header = line[1:].strip()
                                    seq_counter += 1
                                    new_header = ">seq_{}\n".format(seq_counter)
                                    if entry != '':
                                            print(entry)
                                            entry = ''
                                    entry += new_header
                                    renamed_seqs_fh.write("{}\t{}".format(original_header, new_header[1:]))
                            else:
                                    entry += line.strip()
                    print(entry)
            # bc tab_print was True takes this way
            else:
                    if len_print == True:
                            for line in fh:
                                    if line[0] == '>':
                                            if new_header != '':
                                                    seq_len = len(seq)
                                                    print("{0}\t{1}\t{2}".format(new_header, seq_len, seq))
                                                    seq = ''
                                            original_header = line[1:].strip()
                                            seq_counter += 1
                                            new_header = "seq_{}".format(seq_counter)
                                            renamed_seqs_fh.write("{}\t{}\n".format(original_header, new_header))
                                    else:
                                            seq += line.strip()
                            seq_len = len(seq)
                            print("{0}\t{1}\t{2}".format(new_header, seq_len, seq))
                    else:
                            for line in fh:
                                    if line[0] == '>':
                                            if new_header != '':
                                                    print("{0}\t{1}".format(new_header, seq))
                                                    seq = ''
                                            original_header = line[1:].strip()
                                            seq_counter += 1
                                            new_header = "seq_{}".format(seq_counter)
                                            renamed_seqs_fh.write("{}\t{}\n".format(original_header, new_header))
                                    else:
                                            seq += line.strip()
                            print("{0}\t{1}".format(new_header, seq))
    else:
            if tab_print == False:
                    for line in fh:
                            if line[0] == '>':
                                    if entry != '':
                                            print(entry)
                                            entry = ''
                                    entry += line
                            else:
                                    entry += line.strip()
                    print(entry)
            else:
                    if len_print == True:
                            for line in fh:
                                    if line[0] == '>':
                                            if header != '':
                                                    seq_len = len(seq)
                                                    print("{0}\t{1}\t{2}".format(header, seq_len, seq))
                                                    seq = ''
                                            header = line[1:].strip()
                                    else:
                                            seq += line.strip()
                            seq_len = len(seq)
                            print("{0}\t{1}\t{2}".format(header, seq_len, seq))
                    else:
                            for line in fh:
                                    if line[0] == '>':
                                            if header != '':
                                                    print("{0}\t{1}".format(header, seq))
                                                    seq = ''
                                            header = line[1:].strip()
                                    else:
                                            seq += line.strip()
                            print("{0}\t{1}".format(header, seq))
