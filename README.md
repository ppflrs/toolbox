# toolbox
Tools for day to day bioinformatics work


## bamCoverage.py

This script parses (ideally) a BAM file containing alignments to a reference genome. It filters the different alignments by % identity provided with the `-i` parameter.

Some additional info can be obtained by running `bamCoverage.py -h`. 

Because it's an incomplete script it breaks whenever a non-strictly BAM file is supplied, therefore we need to give it a SAM file as input:

`samtools view $BAM_FILE | python bamCoverage.py -i 90 -`

In the command above `bamCoverage.py` takes it's input from the `stdin` (`-`).

The result is written into the file provided as `--outfile` otherwise it will be print in the `stdout` and it's composed by 4 tab-separated columns:

`READ_ID  PCT_ID  BOWTIE2_SCORE ALIGNMENT_POSITION`
