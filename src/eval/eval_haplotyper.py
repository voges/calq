import os.path

basedir = "/data/voges/muenteferi"
datasets = [["ERP001775", "ERR174324.aln_bowtie2.sorted.dupmark.rg.realn.recal"],
            ["NA12878_Garvan_replicate_J", "NA12878_V2.5_Robot_2.aln_bowtie2.sorted.dupmark.rg.realn.recal"],
            ["NA12878-SRX517292", "SRR1238539.aln_bowtie2.sorted.dupmark.rg.realn.recal"]]
subsets = ["3", "11", "20"]

filtersize = ["10", "17", "25"]
filtertype = ["Gauss", "Rectangle"]
quantizerType = ["Uniform", "Lloyd"]

quantSteps = [["2", "8"], ["4", "4"]]

squashed = ["", "--noSquash"]

for dset in datasets:
    for sset in subsets:
        folder = "{}/{}/{}.{}/".format(basedir, dset[0], dset[1], sset)
        file = "{}{}.bam".format(folder, dset[1])
        if os.path.isfile(file):
            print("existing: " + file)
        else:
            print("NOT existing: " + file)



def evaluate (folder):
    pass



