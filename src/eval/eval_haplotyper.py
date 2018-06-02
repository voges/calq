import os.path

basedir = "/data/voges/muenteferi"
datasets = [["ERP001775", "ERR174324.aln_bowtie2.sorted.dupmark.rg.realn.recal"],
            ["NA12878_Garvan_replicate_J", "NA12878_V2.5_Robot_2.aln_bowtie2.sorted.dupmark.rg.realn.recal"],
            ["NA12878-SRX517292", "SRR1238539.aln_bowtie2.sorted.dupmark.rg.realn.recal"]]
subsets = ["3", "11", "20"]

filtersize = ["10", "17", "25"]
filtertype = ["Gauss"]
quantizerType = ["Uniform", "Lloyd"]

quantSteps = [["2", "8"]]

squashed = ["", "--noSquash"]

for dset in datasets:
    for sset in subsets:
        folder = "{}/{}/{}.{}".format(basedir, dset[0], dset[1], sset)
        file = "{}/{}.{}".format(folder, dset[1], sset)
        if not os.path.isfile(file + ".bam"):
            print("NOT existing: " + file + ".bam")
            exit(-1)

        for fsize in filtersize:
            for ftype in filtertype:
                for qtype in quantizerType:
                    for qsteps in quantSteps:
                        for squash in squashed:
                            outfile = "{}.filter{}{}.quant{}_{}{}{}".format(folder,fsize,ftype,qsteps[0], qsteps[1], qtype, squash)
                            print(outfile)



def evaluate (folder):
    pass



