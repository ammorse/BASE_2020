#!/usr/bin/env python3

import argparse, csv
import pandas as pd
import os
import sys
import logging
import logging.config

# McLab Packages
#import mclib_Python as mclib

#  check that design file contains following headers (in order):
#	G1, G2, sampleID, fqName, fqExtension, techRep, readLength
#       fqName is fastq file without extension

def getOptions():
    parser = argparse.ArgumentParser(description='Check design file.  Design file MUST have following columns: G1, G2, sampleID, fqName, fqExtension, techRep, readLength')
    parser.add_argument('-design','--design',dest='design', action='store', required=True, help='Design file containing fq file names and sample ids [Required]')
    parser.add_argument('-design_identifier','--design_identifier',dest='design_identifier', action='store', required=True, help='Design file identifier [Required]')
#    parser.add_argument('-id','--sampleID',dest='sampleID', action='store', required=True, help='Name of the column containing sampleIDs [Required]')
#    parser.add_argument('-fq','--fqName',dest='fqName', action='store', required=True, help='Name of the column containing fastq file names [Required]')
#    parser.add_argument('-g1','--g1',dest='g1', action='store', required=True, help='Name of the column containing G1 names [Required]')
#    parser.add_argument('-g2','--g2',dest='g2', action='store', required=True, help='Name of the column containing G2 names [Required]')
#    parser.add_argument('-rep','--rep',dest='rep', action='store', required=True, help='Name of the column containing tech rep identifiers [Required]')
#    parser.add_argument('-e','--ext',dest='ext', action='store', required=True, help='fastq extension (ex: .fq, .fastq) [Required]')
#    parser.add_argument('-r','--readLen',dest='readLen', action='store', required=True, help='Name of column containing readLength values [Required]')
    # Output data
    parser.add_argument('-d', '--dups', dest='dups', required=False, help='File containing list of duplicate fqNames in design file')
    parser.add_argument('-l', '--logfile', dest='logfile', required=True, help='Name of log file that checks design file')

    args = parser.parse_args()
    return(args)

#check that names of the fastq files in the design file are unique returns duplicate rows if they exist
def fastq_check(design):
    fqName = 'fqName'
    with open(args.dups, 'w') as outfile:
        df = pd.read_csv(design, sep='\t', index_col=None)
      #  print(df)

        if len(df[fqName].unique().tolist()) < len(df[fqName].tolist()):   ##creates lists out of elements in column and compares to a list of unique elements
            dups = df[df.duplicated(fqName, keep=False) == True]           ## creates list of any duplicated rows
            if dups is not None:
                outfile.write('Duplicate check:\tThere are duplicate fqNames in your design file!\t\n---------------------------------\nDuplicated Rows:\n')
                dups.to_csv(outfile, index=False, sep="\t")
                return -1
        else:
            outfile.write('Duplicate check:\tNo duplicate fqNames in your design file\t\n')
            return 0

def columns_check(design, design_identifier):
    columns = ['G1','G2','sampleID','fqName','fqExtension','techRep','readLength']
    with open(args.logfile, 'w') as outfile:
        df = pd.read_csv(design, sep='\t', index_col=None)                
        headers = list(df)
        if headers != columns:
            outfile.write('Headers check:\tERROR: column headers in file ' + design_identifier + ' do not align with order requirements, please check.\t\n---------------------------------\nDetails:' + '\n')
            for col in columns:
                if col not in df:
                    outfile.write('\tError: column header called ' + col + ' does not exist in design file\n')
            return -1
        if headers == columns:
            outfile.write('Headers check:\tColumn headers in file ' + design_identifier + ' align with requirements.\t' + '\n')
            return 0

##run main
def main(args):
    col_match = columns_check(args.design, args.design_identifier)
#    if col_match == -1:
#        raise Exception("ERROR in column headers detected!\n")
    if_dup = fastq_check(args.design)
#    if if_dup == -1:
#        raise Exception("Duplicate rows detected!\n")

if __name__=='__main__':
    args =  getOptions()
    logging.basicConfig(filename=args.logfile,
        filemode='a',
        format='%(asctime)s [%(levelname)s] %(message)s',
        level=logging.DEBUG)
    main(args)


