#!/usr/bin/env python

import sys
import argparse

arg_parser = argparse.ArgumentParser(description='Linearize a FASTA file and print to the standard output.')
arg_parser.add_argument("--tab", "-t", action="store_true", help='Print sequences in tab separated format.')
arg_parser.add_argument("--len", "-l", action="store_true", help='Print lenght of sequence after the sequence header. Only for tab separated format.')
arg_parser.add_argument("input_file",type=str, help='<input_file> full path to fasta file.')

args = arg_parser.parse_args()

infile = args.input_file
tab_print = args.tab
len_print = args.len

if len_print == True:
	if tab_print == False:
		sys.exit("Sorry --len option only valid if used with --tab.")

with open(infile, 'rt') as fh:
	entry = ''
	seq = ''
	header = ''
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
