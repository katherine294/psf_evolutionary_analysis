import pandas as pd

# Load EggNOG annotations
emapper = pd.read_csv("246539E_W163A3B1_chrom.emapper.annotations", sep='\t', comment='#')

# Keep relevant columns
emapper = emapper[['query', 'GOs', 'KEGG_Pathway', 'COG_category']]

# Replace NaN with empty string
emapper = emapper.fillna('')

# Write mapping files
def write_mapping(df, colname, outfile):
    with open(outfile, 'w') as f:
        for _, row in df.iterrows():
            terms = row[colname].replace(',', ';')
            if terms:
                f.write(f"{row['query']}\t{terms}\n")

write_mapping(emapper, 'GOs', 'gene2go.txt')
write_mapping(emapper, 'KEGG_Pathway', 'gene2kegg.txt')
write_mapping(emapper, 'COG_category', 'gene2cog.txt')
