# Nn-Orth_finder (non-orthologue-finder)

"""
NN-Orth_finder (RBH_tools).
Version: alpha
Last modified: 02/11/2020
Author: Cailean Carter

"""

from Bio.KEGG import REST
import pandas as pd  
import time
import requests
import json
import argparse

# parse commandline inputs
parser = argparse.ArgumentParser()

parser.add_argument('rbh', metavar="RBH file path", type=str)

parser.add_argument('annot', metavar="Genome annotation file path", type=str)

parser.add_argument('organism', type=str, 
help="""A KEGG genome identifier for the organism of interest.
KEGG genomes can be found at: https://www.genome.jp/kegg/catalog/org_list.html
If organism is unknown, please provide an identifier for what the organism is suspected to be.
""")

args = parser.parse_args()


try:
    REST.kegg_info(args.organism)
except:
    raise ValueError("KEGG genome (organism) identifier not found. Please provide a valid identifier. See help for info.")


def get_profile(ID):
    """Retrieves the profile page from the 
    KEGG database for a gene or protein"""
    return REST.kegg_get(ID).readlines()

def get_ID(gene, db):
    """Gets the KEGG ID from the KEGG database from 
    which a KEGG profile page can be obtained"""
    result = REST.kegg_find(db, gene).read()

    if "hypothetical protein" in result:
        return False

    ID, *_ = result.split("\t")
    return ID


def get_reaction(profile):
    "Pulls the reaction from the KEGG profile page"

    if "Enzyme" in profile[0] and not "Obsolete" in profile[0]:
        for line in profile[:60]:
            if line.startswith("REACTION"):

                _, *reac = line.strip().split()

                return ' '.join(reac)


def get_pathways(gene, ec, cog):
    "Sends requests to KEGG and COG database for information"

    
    if isinstance(ec, str) and '-' not in ec:
        return get_reaction(
            get_profile(ec)
            )

    elif isinstance(cog, str):

        response = requests.get("http://www.ncbi.nlm.nih.gov/research/cog/api/cogdef/?cog=%s" % cog)
        
        if response.status_code == 200: #make sure request is OK
            
            CogResult = json.loads(response.text)

            if CogResult['results'][0]['pathways']:
                return CogResult['results'][0]['pathways']
                
    
    elif isinstance(gene, str):
        try:
            ID = get_ID(gene, organism)

            if ID:
                profile = get_profile(ID)
                return get_reaction(profile)
        except:
            return

def run(RBH, Annot, organism):

    #store the time at start of script
    start = time.time()

    # read files
    rbh = pd.read_csv(RBH, sep="\t")
    annot = pd.read_csv(Annot, sep="\t")

    # normalise bitscore
    rbh["norm_bitscore"] = rbh.bitscore/rbh.length

    annotIds = list(annot["locus_tag"])
    rbhIds = list(rbh["B_id"])

    # Get locus_tags not present in RBH
    totalexcluded = list(filter(lambda x: x not in rbhIds, annotIds)) 
    # Add locus_tags which were under filter score from RBH
    totalexcluded += list(rbh[rbh["norm_bitscore"] < 1.6]["B_id"]) 


    # Pick out the information from the annotation file for locus_tags
    # which were are non-orthologues. Remove locus_tags that have no
    # database identifiers which can be used to reference to
    # search against the KEGG database or COG database
    nnorths = annot[annot["locus_tag"] \
        .isin(totalexcluded)] \
        .filter(items=["locus_tag", "gene", "EC_number", "COG", "product"]) \
        .dropna(thresh=3)

    # Gets the pathways for all the locus_tags if available from KEGG and COG databases
    nnorths["pathways"] = nnorths.apply(lambda  row: get_pathways(row["gene"], row["EC_number"], row["COG"]), axis=1)


    # Output a table of information from the run.
    statTable = pd.DataFrame.from_dict({

        "Total number of annotations" : annot["locus_tag"].count(),
        "Number of non-orthologues" : nnorths["locus_tag"].count(),
        "Number of pathways found" : nnorths["pathways"].count(),

        "Run time" : "{:.2f} minutes".format( (time.time() - start) /60 )

    }, 
    orient="index")


    print(statTable.to_string(header=False))

    return nnorths

run(args.rbh, args.annot, args.organism)