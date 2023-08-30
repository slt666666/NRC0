import pandas as pd
import numpy as np
import sys
from Bio import SeqIO

specie_class = sys.argv[1]

gene_IDs = []
sequences = []
for record in SeqIO.parse('{}_alignment.fasta'.format(specie_class), 'fasta'):
    gene_IDs.append(record.id)
    sequences.append(list(record.seq))
data = pd.DataFrame(sequences)
data["gene_IDs"] = gene_IDs
data = data.set_index("gene_IDs")

# extract AtZAR1 position
data = data.loc[:, data.loc["AtZAR1", :] != "-"]

AtZAR1_seq = "".join(data.loc["AtZAR1", :].values)
NB_start = AtZAR1_seq.find("NGTDR")
NB_end = AtZAR1_seq.find("LNCRH")

data = data.iloc[:, NB_start:NB_end+5]

# to fasta
with open("{}_alignment_NB_AtZAR1.fasta".format(specie_class), mode='w') as f:
    for i in range(data.shape[0]):
        f.write(">"+data.index[i])
        f.write("\n")
        f.write("".join(data.iloc[i, :].tolist()))
        f.write("\n")

### remove genes which have so many blanks
no_data_count = data == "-"
no_data_count = no_data_count.sum(axis=1)
data = data.loc[no_data_count < 370*0.9, :]

# to fasta
with open("{}_alignment_NB_AtZAR1_rm90.fasta".format(specie_class), mode='w') as f:
    for i in range(data.shape[0]):
        f.write(">"+data.index[i])
        f.write("\n")
        f.write("".join(data.iloc[i, :].tolist()))
        f.write("\n")

### cut based on K of GKTT

with open("{}_xKxx_remove_truncated90.fasta".format(specie_class), mode='w') as f:
    for record in SeqIO.parse("{}_alignment_NB_AtZAR1_rm90.fasta".format(specie_class), 'fasta'):
        id_part = record.id
        desc_part = record.description
        seq = str(record.seq)
        if seq[50] == "K":
            f.write(">"+id_part)
            f.write("\n")
            f.write(seq)
            f.write("\n")

### cut based on G/AxxxxGKT/S of p-loop
with open("{}_GorAxxxxGKTorS_remove_truncated90.fasta".format(specie_class), mode='w') as f:
    for record in SeqIO.parse("{}_alignment_NB_AtZAR1_rm90.fasta".format(specie_class), 'fasta'):
        id_part = record.id
        desc_part = record.description
        seq = str(record.seq)
        if (seq[44] == "G") or (seq[44] == "A"):
            if (seq[49:51] == "GK"):
                if (seq[51] == "T") or (seq[51] == "S"):
                    f.write(">"+id_part)
                    f.write("\n")
                    f.write(seq)
                    f.write("\n")
