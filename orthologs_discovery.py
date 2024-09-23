# In[]
import pandas as pd
import numpy as np

# In[]
# the database stores the gene orthologs, 
# it is downloaded with the link: https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz
gene_ortho_db = pd.read_csv("gene_orthologs", sep = "\t")
# read in the gene information, the information and readme under the folder: https://ftp.ncbi.nlm.nih.gov/gene/DATA/
gene_info = pd.read_csv("gene_info", sep = "\t") # including the readme information
# gene2ensemble = pd.read_csv("gene2ensembl", sep = "\t")

# column ``relationship'' has only one value ``Ortholog''
gene_ortho_db = gene_ortho_db.drop(["relationship"], axis = 1)
# each row store one pair of orthologs, tax_id is the species id, and gene_id is the gene id.
gene_ortho_db.columns = ["tax_id_1", "GeneID_1", "tax_id_2", "GeneID_2"]
# full name of tax id is taxonomy id
# species_taxid = pd.DataFrame(columns = ["species", "tax_id"])
# species_taxid["species"] = ["human", "marmoset", "mouse", "armadillo_9_banded", "opossum", "naked_mole_rat", "squirrel_arctic_ground", "tree_shrew", "pig_tailed_macaque", "rat", "pig", "ferret", "rabbit", "dog", "cat", "rhesus_macaque", "african_green_monkey", "chimpanzee", "gorilla", "galago", "olive_baboon", "squirrel_monkey", "owl_monkey", "mouse_lemur"]
# species_taxid["tax_id"] = [9606, 9483, 10090, 9361, 13616, 10181, 9999, 37347, 9545, 10116, 9823, 9669, 9986, 9615, 9685, 9544, 60711, 9598, 9595, 30611, 9555, 39432, 37293, 30608]
# species_taxid.to_csv("species_taxid.txt", sep = "\t")
species_taxid = pd.read_csv("species_taxid.txt", sep = "\t", index_col = 0)

def ortho_species(species1_tax_id, species2_tax_id, ortho_db, gene_info):
    """\
    Description:
    --------------
        Subsetting the orthologs database to find the orthogs terms only relevant to species 1 and species 2, and return cross-species orthologs
    
    Parameters:
    -------------
        species1_tax_id: the taxonomy id of species 1
        species2_tax_id: the taxonomy id of species 2
        ortho_df: the gene orthologs dataframe, including columns ``tax_id_1'', ``tax_id_2'', ``GeneID_1'' and ``GeneID_2''
        gene_info: the information of genes
    Returns:
    ------------
        ortho_db_f: the selected orthologs dataframe including only species 1 and species 2
        ortho_mtx: the binary orthologs matrix (including one-one, one-many) between species 1 (row) and species 2 (column), the index and column names are gene ids.
    """
    ortho_db_sp12 = ortho_db[(ortho_db["tax_id_1"] == species1_tax_id) & (ortho_db["tax_id_2"] == species2_tax_id)]
    ortho_db_sp21 = ortho_db[(ortho_db["tax_id_1"] == species2_tax_id) & (ortho_db["tax_id_2"] == species1_tax_id)]
    # make sure the first species is always species 1
    ortho_db_sp21 = ortho_db_sp21.iloc[:, [2,3,0,1]]
    ortho_db_sp21.columns = ortho_db_sp12.columns.values
    ortho_db_f = pd.concat([ortho_db_sp12, ortho_db_sp21], axis = 0, ignore_index = True)

    # transform the gene id information into the gene symbol
    gene_info_sp1 = gene_info[gene_info["#tax_id"] == species1_tax_id]
    gene_info_sp2 = gene_info[gene_info["#tax_id"] == species2_tax_id]
    gene_info_sp1.index = gene_info_sp1["GeneID"].values
    gene_info_sp2.index = gene_info_sp2["GeneID"].values
    
    ortho_db_f["Symbol_1"] = gene_info_sp1.loc[ortho_db_f["GeneID_1"].values, "Symbol"].values
    ortho_db_f["Symbol_2"] = gene_info_sp2.loc[ortho_db_f["GeneID_2"].values, "Symbol"].values

    # build a binary orthologs matrix for each pair of species, pros: store one-to-one and one-to-many
    ortho_db_tmp = ortho_db_f[["Symbol_1", "Symbol_2"]]
    ortho_db_tmp["value"] = 1
    ortho_mtx = ortho_db_tmp.pivot(index = "Symbol_1", columns = "Symbol_2", values = "value").fillna(0).astype(int)    
    ortho_mtx.columns = ortho_mtx.columns.values
    ortho_mtx.index = ortho_mtx.index.values    
    ortho_mtx = ortho_mtx.astype(pd.SparseDtype(int, fill_value=0))
    return ortho_db_f, ortho_mtx

# In[]
human_taxid = species_taxid.loc[species_taxid["species"] == "human", "tax_id"].values[0]
mouse_taxid = species_taxid.loc[species_taxid["species"] == "mouse", "tax_id"].values[0]
marmoset_taxid = species_taxid.loc[species_taxid["species"] == "marmoset", "tax_id"].values[0]
ortho_human_mouse_df, ortho_human_mouse_mtx = ortho_species(human_taxid, mouse_taxid, gene_ortho_db, gene_info)
ortho_human_marmoset_df, ortho_human_marmoset_mtx = ortho_species(human_taxid, marmoset_taxid, gene_ortho_db, gene_info)

# ------------------------------------------------------------------------------------------------------------
# NOTE: the gene gtf information does not match the ortho above
# download the ncbi gtf file for each species of interest.
# On the ncbi website, the gtf (gene annotation) file can be downloaded with the link https://ftp.ncbi.nlm.nih.gov/genomes/refseq
# For example, to download the marmoset (formal name Callithrix_jacchus) gtf with:
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Callithrix_jacchus/latest_assembly_versions/GCF_011100555.1_mCalJa1.2.pat.X/GCF_011100555.1_mCalJa1.2.pat.X_genomic.gtf.gz
# for mouse:
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz  
# for human:
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
#   
# In[]

# %%
