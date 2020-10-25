#! /usr/bin/Rscript
import subprocess
import os
import pandas as pd
import argparse
import re

def getOptions():
    parser = argparse.ArgumentParser(description= "Run bayesian model")
    parser.add_argument("-datafile","--datafile",dest="datafile",action="store",required=True,help="Provide path for input datafile [CSV]")
    parser.add_argument("-compID","--compID",dest="compID",action="store",required=True,help="Provide comparison identifier")
    parser.add_argument("-comparate_1","--comparate_1",dest="comparate_1",action="store",required=True,help="Provide comparate 1")
    parser.add_argument("-comparate_2","--comparate_2",dest="comparate_2",action="store",required=True, help="Provide comparate 2")
    parser.add_argument("-datafile2","--datafile2",dest="datafile2",action="store",required=True,help="Provide temp path for created datafile with headers for Bayesian")
    parser.add_argument("-cond","--cond", dest="cond", action="store", required=False, help="Number of conditions")
    parser.add_argument("-workdir","--workdir", dest="workdir", action="store", required=True, help="Path to R code")
    parser.add_argument("-routput","--routput", dest="routput", action="store", required=False, help="Optional R output file")
    parser.add_argument("-subpath","--subpath", dest="subpath", action="store", required=True, help="Bayesian R script subprocess path")
    parser.add_argument("-iterations","--iterations", dest="iterations", action="store", required=False, help="Optional number of iterations (default 100000)")
    parser.add_argument("-warmup","--warmup", dest="warmup", action="store", required=False, help="Optional number of warmup (default 10000)")
    parser.add_argument("-o","--output",dest="output",action="store", required=True,help="Output path for merged file[CSV]")
    args=parser.parse_args()
    return(args)
    args = parser.parse_args()

def main():
    args = getOptions()
    
    ##(1) Parsing datafile to extract rows with sampleID specified in design file, set c1 and c2

    ##Standardized Paths##
    args.output = os.path.abspath(args.output)
    args.routput = os.path.abspath(args.routput)

    #Make variable for number of conditions
    compnum=args.cond

    #Make variable for working directory (stan and wrapper in same place)
    workdir = args.workdir

    comparison=args.compID
    c1=args.comparate_1
    c2=args.comparate_2

#    print(comparison)
#    print(c1)
#    print(c2)

#    row_list = []
#    row_list = [c1, c2, comparison]

#    myDict = {}
#    myDict = 
#    for index, sample in 
#    comparison = comparison_list[key]
#    c1= c1_list[key]
#    c2= c2_list[key]

    infileName = "bayesian_input_" + comparison + ".csv"

    infile=pd.read_csv(os.path.join(args.datafile, infileName))
    infile.set_index('FEATURE_ID')


    ## add comparison column (last) 
    infile['comparison']=comparison
    
    datafile2 = args.datafile2 + comparison + '_temp.csv'

    infile.to_csv(datafile2, na_rep = 'NA', index=False)

    rout = comparison + '_r_out.csv'

    routput=os.path.join(args.routput, rout)

    ## (2) Calls subprocess to run R script where args1 is the input csv and args2 is the output path for NBmodel.R##
    rscript = subprocess.call(['Rscript', args.subpath, datafile2, routput, compnum, workdir, args.iterations, args.warmup])

    ## (3) Format input from Rscript and get list of default header names, change headers back to actual comparates from c1 and c2
    df2 = pd.read_csv(routput)

    print("Printing datafile2 - input into R code")
    print(datafile2)

    headers_out =list(df2.columns.values)

    for a in range(len(headers_out)):
        if 'c1' in headers_out[a]:
            headers_out[a] = headers_out[a].replace('c1', c1)
        if 'c2' in headers_out[a]:
            headers_out[a] = headers_out[a].replace('c2', c2)

    print('printing new headers_out')
    print(headers_out)

    df2.columns = headers_out
    
    ##Write to new CSV##
    outfile = 'bayesian_out_' + comparison + '.csv'
    output=os.path.join(args.output, outfile)

    df2.to_csv(output, na_rep = 'NA', index=False)

if __name__=='__main__':
    main()


