import os.path

basedir = "/data/voges/muenteferi"
datasets = [["ERP001775", "ERR174324.aln_bowtie2.sorted.dupmark.rg.realn.recal"],
            ["NA12878_Garvan_replicate_J", "NA12878_V2.5_Robot_2.aln_bowtie2.sorted.dupmark.rg.realn.recal"],
            ["NA12878-SRX517292", "SRR1238539.aln_bowtie2.sorted.dupmark.rg.realn.recal"]]
subsets = ["3"]

filtersize = ["10"]
filtertype = ["Gauss"]
quantizerType = ["Uniform"]

quantSteps = [["2", "8"]]

squashed = [""]

install_path = "/project/dna/install"
samtools = install_path + "/samtools-1.3/bin/samtools"
calqPath = "/home/muenteferi/Dokumente/calqBuild/calq"
referencePath = "/data/voges/muenteferi/GATK_bundle-2.8-b37/human_g1k_v37.20.fasta"

for dset in datasets:
    for sset in subsets:
        folder = "{}/{}/{}.{}".format(basedir, dset[0], dset[1], sset)
        filename = "{}.{}".format(dset[1], sset)
        filepath = "{}/{}".format(folder, filename)
        if not os.path.isfile(filepath + ".bam"):
            print("NOT existing: " + filepath + ".bam")
            exit(-1)

        if not os.path.isfile(filepath + ".sam"):
            print("Converting to sam: " + file + ".bam" )
            command = "{} view -h -o {}.sam {}.bam".format(samtools, filepath, filepath)
            os.system(command)

        for fsize in filtersize:
            for ftype in filtertype:
                for qtype in quantizerType:
                    for qsteps in quantSteps:
                        for squash in squashed:
                            outfolder = "{}.calq-haplo.filter{}{}.quant{}_{}{}{}".format(folder,fsize,ftype,qsteps[0], qsteps[1], qtype, squash)
                            outfile = "{}.calq-haplo.filter{}{}.quant{}_{}{}{}".format(outfolder + "/" + filename, fsize, ftype, qsteps[0], qsteps[1], qtype, squash)
                            if not os.path.isdir(outfolder):
                                print("Creating dir: " + outfolder)
                               # os.system("mkdir " + outfolder)
                            calqCommand = "{} -q Illumina-1.8+ -p 2 -b 10000 {}  -o {} -r {] --quantizerType {} " \
                                          "--filterType {} --quantizationMin {} --quantizationMax {} --filterSize {} {}".\
                                format(calqPath, filepath + ".sam", outfile, referencePath, qtype, ftype, qsteps[0],
                                       qsteps[1], fsize, squash)
                            print(calqCommand)







