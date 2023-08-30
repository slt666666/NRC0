# import library
import pandas as pd
import numpy as np
from Bio import SeqIO
import glob

# make gff format data for NLRs
species = ["Coffea_canephora", "Daucus_carota", "Solanum_lycopersicum", "Ipomoea_trifida"]
gff_type = ["gene", "gene", "mRNA", "mRNA"]
gff_files = glob.glob("*gff*") # read path to gff files

for i, specie in enumerate(species):
    gff = pd.read_csv(gff_files[i], comment="#", sep="\t", header=None)
    gff = gff[gff.iloc[:, 2] == gff_type[i]]
    NLR_list = [record.id for record in SeqIO.parse("{}_NLR.fasta".format(specie), 'fasta')] # read NLR pep ids from fasta files
    NLR_gff = pd.DataFrame()
    for each_NLR in sorted(NLR_list):
        NLR_gff = pd.concat([NLR_gff, gff[gff.iloc[:, 8].str.contains(each_NLR)]])
    NLR_gff.columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    NLR_gff["pep_id"] = sorted(NLR_list)
    NLR_gff.to_csv("{}_NLR.csv".format(specie), index=None) # save gff information as csv format

# identify clustered pair NLRs
def make_clustered_pair_NLRs(specie, threshold_distance):
    NLR_gff = pd.read_csv("{}_NLR.csv".format(specie))
    clustered_pair_NLRs = []
    for chrom in NLR_gff.seqname.unique():
        chrom_gff = NLR_gff[NLR_gff["seqname"] == chrom]
        for each_id in chrom_gff.pep_id.values:
            each_id_gff = chrom_gff[chrom_gff["attribute"].str.contains(each_id)]
            clustered_gff = chrom_gff[(abs(chrom_gff["start"] - each_id_gff.end.values[0]) <= threshold_distance) | (abs(each_id_gff.start.values[0] - chrom_gff["end"]) <= threshold_distance)]
            if clustered_gff.shape[0] > 1:
                clustered_gff = clustered_gff[clustered_gff["pep_id"] != each_id]
                for clustered_NLR in clustered_gff.pep_id.values:
                    distance = min(abs(clustered_gff[clustered_gff["pep_id"] == clustered_NLR].start.values[0] - each_id_gff.end.values[0]), abs(each_id_gff.start.values[0] - clustered_gff[clustered_gff["pep_id"] == clustered_NLR].end.values[0]))                      
                    clustered_pair_NLR = sorted([each_id, clustered_NLR])
                    clustered_pair_NLR.append(distance)
                    clustered_pair_NLRs.append(clustered_pair_NLR)

    clustered_pair_NLRs = np.unique(clustered_pair_NLRs, axis=0)
    return clustered_pair_NLRs

threshold_distance = 50000
clustered_pair_NLRs = make_clustered_pair_NLRs("Solanum_lycopersicum", threshold_distance)
print(clustered_pair_NLRs)