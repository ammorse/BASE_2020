#!/usr/bin/env python
import csv, sys, argparse, os, itertools, operator, collections
import pandas as pd
import numpy as np
from csv import DictReader
csv.field_size_limit(1000000000)

### check read numbers into and out of sam compare
### outputs a file 


# Parse command line arguments
parser = argparse.ArgumentParser(description='check read numbers into and out of sam compare, must be within minimum unique reads and sum of uniq reads from both summary files')
parser.add_argument('-s1','--summary1',dest='summary1', action='store', required=True, help='The sam summary file containing read counts after dropping [Required]')
parser.add_argument('-s2','--summary2',dest='summary2', action='store', required=True, help='The sam summary file containing read counts after dropping [Required]')
parser.add_argument('-ase_names','--ase_names',dest='ase_names', action='store', required=True, help='fastq filename [Required]')
parser.add_argument('-a','--ase',dest='ase',action='store',required=True, help='The ase totals file containing read counts generated from sam compare script [Required]')
parser.add_argument('-o','--out', dest='out', action='store', required=True, help='Output file containing check info [Required]')
args = parser.parse_args()
  

       

B1 = []
B2 = []
### Open bwa
with open(args.summary1,'r') as bwa_table:
    B1 = csv.reader(bwa_table, delimiter= '\t')
    next(B1)
    for row in B1:
        uniq_b1 = int(row[1])
        #total reads after dropping

with open(args.summary2, 'r') as bwa2_table:
    B2 = csv.reader(bwa2_table, delimiter= '\t')
    next(B2)
    for row in B2:
        uniq_b2 = int(row[1])
        #total reads after dropping

## read # in sam file should be between greater of uniq_b1 or uniq_b2 - this should be the same # as in the ase_totals table
sumReads=uniq_b1 + uniq_b2
minReads=min(uniq_b1, uniq_b2)

print(sumReads)
print(minReads)

with open(args.ase, 'r') as ase_table:
    df = pd.read_csv(ase_table, sep='\t')

    count_tot=df['Count totals:'].iloc[len(df)-1]
    print(count_tot)

if int(minReads) <= int(count_tot) <= int(sumReads):
    flag_readnum_in_range = 1
else:
    flag_readnum_in_range = 0


name = os.path.basename(args.ase_names)
#bname = os.path.basename(args.fq)
#name  = os.path.splitext(bname)[0]

## counts in ase file should be betweeen minReads and the sum of uniq mapped reads 
## open file to write to
with open(args.out, 'w') as outfile:
    spamwriter=csv.writer(outfile, delimiter='\t')
    first_row = True
    
    if first_row: 
        spamwriter.writerow(['fqName', 'min_uniq_g1_uniq_g2', 'sum_uniq_g1_uniq_g2', 'total_counts_ase_table', 'flag_readnum_in_range'])
        first_row = False

        row_items = [name, minReads, sumReads, count_tot,flag_readnum_in_range]
    spamwriter.writerow(row_items)

