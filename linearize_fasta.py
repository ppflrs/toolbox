#!/usr/bin/env python

import sys
import argparse
import os

arg_parser = argparse.ArgumentParser(description='Linearize a FASTA file and print to the standard output.')
arg_parser.add_argument("--tab", "-t", action="store_true", help='Print sequences in tab separated format.')
arg_parser.add_argument("--len", "-l", action="store_true", help='Print lenght of sequence after the sequence header. Only for tab separated format.')
arg_parser.add_argument("--min_len", "-ml", type=int, help='Prints sequences wich lenght is bigger than --min-len.')
arg_parser.add_argument("--rename", "-r", type=str, help='<tsv_file> Rename the headers to "seq_1, seq_2, etc." format. The original headers are stored in tsv_file.')
arg_parser.add_argument("input_file",type=str, help='<input_file> full path to fasta file.')

args = arg_parser.parse_args()

infile = args.input_file
tab_print = args.tab
len_print = args.len
rename_seqs_file = args.rename
min_len = args.min_len

if rename_seqs_file != None:
    rename_seqs = True
    if rename_seqs_file == infile:
        sys.exit("Warning: This will delete the input file, please change the renamed headers mapping file.")
    if os.path.exists(rename_seqs_file):
        sys.exit("Error: The file to store the renamed headers already exists.")
    else:
        open(rename_seqs_file, 'a')
else:
    rename_seqs = False
    rename_seqs_file = False
seq_counter = 0

if min_len == None:
    # min_len = False
    min_len = 0
elif min_len < 0:
    sys.exit("Exiting: --min_len should be a positive integer.")

if len_print == True:
    if tab_print == False:
        sys.exit("Sorry --len option only valid if used with --tab.")

def print_entry(entry, tab_print, len_print, min_len, seq_counter, rename_seqs_file):
    header = entry["header"]
    seq = entry["seq"]
    seq_len = len(seq)

    if rename_seqs_file:
        with open(rename_seqs_file, 'at') as fh:
            fh.write("{}\tseq_{}\n".format(header, seq_counter))
            header = "seq_{}".format(seq_counter)
    if seq_len >= min_len:
        if tab_print:
            if len_print:
                print("\t".join(map(str, [header, seq_len, seq])))
            else:
                print("\t".join(map(str, [header, seq])))
        else: # no tab_print
            print(">{}\n{}\n".format(header, seq))
    else: # seq_len < min_len
        return # don't print anything

def parse_file(fh, tab_print, len_print, min_len, rename_seqs_file, seq_counter):
    entry = {"header":"", "seq":""}
    header = entry["header"]
    for line in fh:
            if line[0] == '>':
                    if entry["header"] != '':
                            print_entry(entry, tab_print, len_print, min_len, seq_counter, rename_seqs_file)
                            entry = {"header":"", "seq":""}
                    entry["header"] = line[1:].strip()
                    seq_counter += 1
            else:
                    entry["seq"] += line.strip()
    print_entry(entry, tab_print, len_print, min_len, seq_counter, rename_seqs_file)
    return

if infile == '-':
    fh = sys.stdin.readlines()
    parse_file(fh, tab_print, len_print, min_len, rename_seqs_file, seq_counter)
else:
    with open(infile, 'rt') as fh:
        parse_file(fh, tab_print, len_print, min_len, rename_seqs_file, seq_counter)

