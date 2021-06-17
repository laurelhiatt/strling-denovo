import argparse
import pandas as pd
import numpy as np
import peddy ### must be in python version 3.7 for peddy to actually work

def has_parents(sample):
    """Check if Peddy sample has both parents in ped file"""
    if sample.mom is not None and sample.dad is not None:
        return True
    return False
### this step is important so that we can identify the trios that we are working with

def get_args():
    """Incorporating argparse into the code for interchangeable arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--outliers",
        help="input outlier name")
    parser.add_argument("--ped",
        help="input ped file")
	### out will just be the name of my output file... turn up
    parser.add_argument("--out",
        help="outputfile")
    return parser.parse_args()

def expandorama(df,kid,mom,dad, writeHeader = True):
    """Generate .tsv file(s) with pedigree input and STRling data that output max allele differences and if
	kid has expansion greater than both mom and dad"""
    args = get_args()
    df["allelecomp"] = df[["allele1_est", "allele2_est"]].max(axis=1)
    df['allelecomp'] = df['allelecomp'].replace(np.nan, 0)
    ###first, we want to get the max allele value and replace any NaNs with 0 so we can math correctly

    dfkid = df.loc[df['sample'] == kid]
    dfmom = df.loc[df['sample'] == mom]
    dfdad = df.loc[df['sample'] == dad]
	### this is how we match our pedigree samples to our data frame samples, with the sample IDs

    dfkid = dfkid.rename(columns={"allelecomp": "allele_kid"})
    dfdad = dfdad.rename(columns={"allelecomp": "allele_dad"})
    dfmom = dfmom.rename(columns={"allelecomp": "allele_mom"})
    ### since we know that all of the alleles are composite, we rename them to tell apart the trio members

    drop_from_dkid= (['allele1_est', 'allele2_est','spanning_reads', 'spanning_pairs', 'left_clips', 'right_clips', 'unplaced_pairs', 'sum_str_counts', 'sum_str_log', 'depth', 'outlier'])
    drop_from_parents = (['left', 'right', 'chrom', 'sample', 'p', 'p_adj', 'repeatunit']) + drop_from_dkid
    dfkid = dfkid.drop(drop_from_dkid, axis=1)
    dfmom = dfmom.drop(drop_from_parents, axis=1)
    dfdad = dfdad.drop(drop_from_parents, axis=1)
    ### we are dropping as many columns as we can for a clean output, while still getting essential information

    kiddad = dfkid.merge(dfdad, on= 'locus')
    kiddadmom = kiddad.merge(dfmom, on= 'locus')
	### we are merging the dataframes by locus so we can easily subtract the columns, and have a clean output by locus

    kiddadmom = kiddadmom.assign(kiddeldad=kiddadmom['allele_kid'] - kiddadmom['allele_dad'])
    kiddadmom = kiddadmom.assign(kiddelmom=kiddadmom['allele_kid'] - kiddadmom['allele_mom'])
	###we are creating a new column that is the difference between child and parent, which gives an idea of the expansions

    kiddadmom['novel_amp'] = (kiddadmom['allele_kid']-kiddadmom['allele_dad']>0) & (kiddadmom['allele_kid']-kiddadmom['allele_mom']>0)
	### we make a new column where the difference between child and parent is positive for both, prints True; these are candidate expansions

    novel_amp_reads = kiddadmom.novel_amp.value_counts()
    ### We get a true/false count per trio. Neat!
    my_small_df = (kiddadmom, novel_amp_reads, 'Kid, mom, and dad sample IDs are', kid, mom, dad)
	### I just kinda like this dataframe, not super useful but yeah
    if writeHeader is True:
        kiddadmom.to_csv(args.out, mode='a', sep='\t', header=True, index=False)
        writeHeader = False
    else:
        kiddadmom.to_csv(args.out, mode='a',sep='\t', header=False, index=False)
    ###kiddadmom.to_csv(args.out, mode = 'a', sep='\t', header = write_header, index = False)
    print(novel_amp_reads, kid) ### to summarize expansions and list the child of trio per dataset
    return my_small_df ### if I want the dataframe as an object, although it is saved to the composite file

def main():    ###match below or else
    args = get_args()
    df = pd.read_table(args.outliers, delim_whitespace = True, dtype = {'sample' : str}, index_col = False)
    ped = peddy.Ped(args.ped, 'Paternal_ID' == str) ### import the ped file through a peddy function
    ###this is where we input our STRLing outlier data, super exciting!
    with open(args.out, 'w') as newfile:
            pass
    writeHeader = True
    for sample in ped.samples():
            if has_parents(sample):
                expandorama(df, sample.sample_id, sample.maternal_id,sample.paternal_id, writeHeader)
                writeHeader = False ###don't want to keep writing header
    ###kiddadmom.to_csv(args.out, header=kiddadmom.colnames, sep='\t' index=False)

   ### Theoretically at the end, I will have something like the below.
###python strling-denovo.py --outliers STR.tsv --ped families.ped --out my_output.tsv


if __name__ == "__main__": ### don't change this word
	main()  ### this line is a free for all, but by convention we write a function called main. entry point.
