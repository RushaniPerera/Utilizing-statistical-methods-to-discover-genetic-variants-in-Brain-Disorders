import argparse,gzip,os
import pandas as pd
import numpy as np

pd.options.mode.chained_assignment = None # Getting rid of a false positive warning

parser = argparse.ArgumentParser(description='This is a script that does SNP to gene mapping. First argument is the SNP file and second is the gene location file')

# Add arguments
parser.add_argument('snpfile', type=str, help='path to the SNP file')
parser.add_argument('genefile', type=str, help='path to the gene location file')
parser.add_argument('--sep', type=str, default='\s+', help='separator for the gene file, default is spaces')
parser.add_argument('--pval', type=float, default="1", help='optional critical p-value for SNPs, defaults to 1.0')
parser.add_argument('--cs2gdir', type=str, default="", help='optional directory for cs2g data. \
                    If not given then output will only contain SNPs that can be directly mapped.')
parser.add_argument('--cs2g_cutoff', type=float, default="0.0", help='optional cutoff for cs2g data. \
                    If not given then default is 0')

# Parse the arguments
args = parser.parse_args()


genelistfile = args.genefile
genelistfile_sep = args.sep
top_file = args.snpfile
critical_p = args.pval
cs2gdir = args.cs2gdir
cs2g_cutoff = args.cs2g_cutoff

def break_on_chromosome(df):
    return [pd.DataFrame( df.loc[ df[df.columns[0]]==chromosome ,:] ) for chromosome in pd.unique( df[df.columns[0]] )]

def cS2G_scored_SNP_data(cs2gdir):
    # Initializing variables
    chromsnps = [[] for i in range(22)]
    cs2gdirfiles = os.listdir(cs2gdir)
    snpfiles = [i for i in cs2gdirfiles if 'cS2G' in i]

    # Get the SNP locations
    with gzip.open(os.path.join(cs2gdir,"allsnps.txt.gz"),'rt') as file:
        all_snp_locations = pd.read_csv(file,sep=' ',header=None)
    all_snp_locations.columns = ['SNP','CHR','LOC']

    # Read in SNP scores and merge with relevent SNP locations
    for file in snpfiles:

        if file.split('.')[-1] == 'gz':
            with gzip.open(os.path.join(cs2gdir,file),'rt') as filebuffer:
                datafile = pd.read_csv(filebuffer,sep='\t',header=0)
            datafile = datafile.merge(all_snp_locations,on='SNP')
            chromsnps[int(file.split('.')[1])-1] = datafile

    return pd.concat(chromsnps)

def match_snp_location(gene_locs, snp_locs, critical_p=1,cs2g_data=None,cs2g_cutoff=0):
    """
    This function takes in a pandas dataframe of gene locations and SNP locations and determines which gene each SNP can be mapped to.
    This assumes that both dataframes relate to the ***same chromosome***. the 3rd and 4th columns of the gene locations dataframe
    should be the start and end index of the gene

    
    """

    starts = gene_locs[gene_locs.columns[2]].to_numpy()
    ends = gene_locs[gene_locs.columns[3]].to_numpy()
    geneids = []
    genenames = []
    snp_locs = snp_locs.loc[:,['chr','location','pvalue']]

    if cs2g_data is not None:
        # Check which SNPs in cs2g data are in the top file
        snpsintop = np.isin(snp_locs['location'].values,cs2g_data['LOC'].values)
        snp_cs2g = snp_locs.loc[snpsintop,['chr','location','pvalue']]
        snp_cs2g = snp_cs2g.merge(cs2g_data,left_on='location',right_on='LOC',how='left') # Left join to add columns from cs2g file to identified SNPs
        snp_cs2g = snp_cs2g.loc[snp_cs2g['cS2G'] > cs2g_cutoff,['chr','location','pvalue','GENE','SNP','cS2G']]
        snp_locs = snp_locs.iloc[~snpsintop,]

    # iterate row by row
    lines = len(snp_locs)
    for _,snp in snp_locs.iterrows():
        if float(snp['pvalue']) > critical_p: # check if snp p-value is smaller than critical
            continue
        snplocation = int(snp['location'])
        geneloc = np.argmax(starts > snplocation) -1

        if ends[geneloc] < snplocation:
            geneids.append("NA")
            genenames.append("NA")
        else:
            geneids.append(gene_locs.iloc[geneloc,4])
            genenames.append(gene_locs.iloc[geneloc,5])
    
    snp_locs['gene_ID'] = geneids
    snp_locs['gene_name'] = genenames

    if cs2g_data is not None:
        snp_cs2g = snp_cs2g.rename(columns={'GENE':'gene_name'})
        snp_locs = pd.concat([snp_locs, snp_cs2g], ignore_index=True)
    return snp_locs

if __name__=="__main__":    
    geneframe = pd.read_csv(genelistfile, sep=genelistfile_sep, header=None)
    chrframes = break_on_chromosome(geneframe)
    
    if cs2gdir != '':
        cs2g_data = cS2G_scored_SNP_data(cs2gdir)

    if len(top_file.split('\n')) > 1:

        for topfile in top_file.split('\n'):
            output_path = os.path.split(topfile)[0]
            with open(topfile) as file:

                header = file.readline() # Skip header
                columns = file.readline()[1:].strip().split(',') # Get column names

                snpfile = pd.read_csv(file)
                snpfile.columns = columns
            
            snpframe = break_on_chromosome(snpfile)

            # Make the new SNP dataframe with geneIDs added
            if cs2gdir != '':
                newsnps = [match_snp_location(chrom,snps,critical_p,cs2g_data, cs2g_cutoff) for chrom,snps in zip(chrframes,snpframe)]
                newsnps = pd.concat(newsnps)
            else:
                newsnps = [match_snp_location(chrom,snps,critical_p) for chrom,snps in zip(chrframes,snpframe)]
                newsnps = pd.concat(newsnps)

            disease = os.path.split(output_path)[-1]
            if cs2gdir != '':
                disease = 'cS2G_' + disease

            with open(os.path.join(output_path,disease + "_geneIDs.txt"),'w') as file:
                print(f"outputting to {output_path}")
                file.write(','.join(newsnps.columns) + '\n')
                file.write(newsnps.to_csv(index=False,header=False,lineterminator='\n'))
            
            print(f"done with file: {topfile}")
