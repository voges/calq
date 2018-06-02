basedir = "/data/voges/muenteferi"
datasets = ["ERP001775/ERR174324.aln_bowtie2.sorted.dupmark.rg.realn.recal.",
            "NA12878_Garvan_replicate_J/NA12878_V2.5_Robot_2.aln_bowtie2.sorted.dupmark.rg.realn.recal.",
            "NA12878-SRX517292/SRR1238539.aln_bowtie2.sorted.dupmark.rg.realn.recal."]
subsets = ["3", "11", "20"]

for dset in datasets:
    for sset in subsets:
        file = "{}/{}.{}/".format(basedir, dset, sset)
        print(file)


