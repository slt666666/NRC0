mafft --anysymbol NLR.fasta > NLR_alignment.fasta

python3 modules/edit_alignment_for_phylogeny.py NLR
