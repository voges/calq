import os.path

basedir = "/data/voges/muenteferi"
#datasets = [["ERP001775", "ERR174324.aln_bowtie2.sorted.dupmark.rg.realn.recal"],
#            ["NA12878_Garvan_replicate_J", "NA12878_V2.5_Robot_2.aln_bowtie2.sorted.dupmark.rg.realn.recal"],
#            ["NA12878-SRX517292", "SRR1238539.aln_bowtie2.sorted.dupmark.rg.realn.recal"]]
datasets = [["ERP001775", "ERR174324.aln_bowtie2.sorted.dupmark.rg.realn.recal"]]
subsets = ["20"]

filtersize = ["10"]
filtertype = ["Gauss"]
quantizerType = ["Uniform"]

quantSteps = [["2", "8"]]

squashed = [""]

install_path = "/project/dna/install"
samtools = install_path + "/samtools-1.3/bin/samtools"
calqPath = "/home/muenteferi/Dokumente/calqBuild/calq"
referencePath = "/data/voges/muenteferi/GATK_bundle-2.8-b37/human_g1k_v37.20.fasta"
replacePath = "/home/muenteferi/Dokumente/calq/src/ngstools/replace_qual_sam.py"
platypusPath = "/home/muenteferi/Dokumente/calq/src/variant_calling_pipelines/Platypus.sh"
GATK_HF_Path = "/home/muenteferi/Dokumente/calq/src/variant_calling_pipelines/GATK_HF.sh"
GATK_VQSR_Path = "/home/muenteferi/Dokumente/calq/src/variant_calling_pipelines/GATK_VQSR.sh"
HAPPY_prefix = "/home/muenteferi/Dokumente/calq/src/variant_calling_pipelines/do_happy+reppy_chr"

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
                            # Create folders and SAM input
                            outfolder = "{}.calq-haplo.filter{}{}.quant{}_{}{}{}".format(folder,fsize,ftype,qsteps[0], qsteps[1], qtype, squash)
                            outfile = "{}.calq-haplo.filter{}{}.quant{}_{}{}{}.out".format(outfolder + "/" + filename, fsize, ftype, qsteps[0], qsteps[1], qtype, squash)
                            if not os.path.isdir(outfolder):
                                print("Creating dir: " + outfolder)
                                os.system("mkdir " + outfolder)

                            # Compress
                            calqCommand = "{} -q Illumina-1.8+ -p 2 -b 10000 {}.sam  -o {} -r {} --quantizerType {} " \
                                          "--filterType {} --quantizationMin {} --quantizationMax {} --filterSize {} {}".\
                                format(calqPath, filepath, outfile, referencePath, qtype, ftype, qsteps[0],
                                       qsteps[1], fsize, squash)
                            print(calqCommand + "\n")
                            os.system(calqCommand)

                            # Decompress
                            calqDeCommand = "{} -d -s {}.sam {} -o {}.qual".format(calqPath, filepath, outfile, outfile)
                            print(calqDeCommand + "\n")
                            os.system(calqDeCommand)

                            # Insert quality values
                            replaceCommand = "python {} {}.sam {}.qual 1> {}.sam".format(replacePath, filepath, outfile, outfile)
                            print (replaceCommand + "\n")
                            os.system(replaceCommand)

                            # Create bam
                            convertCommand = "{} view -bh {}.sam > {}.bam".format(samtools, outfile, outfile)
                            print(convertCommand + "\n")
                            os.system(convertCommand)

                            # Index bam
                            indexCommand = "{} index -b {}.bam {}.bai".format(samtools, outfile, outfile)
                            print(indexCommand + "\n")
                            os.system(indexCommand)

                            # Platypus
                            PlatypusCommand = "{} 12 {}.bam {}".format(platypusPath, outfile, sset)
                            print(PlatypusCommand + "\n")
                            os.system(PlatypusCommand)

                            # GATK_HF
                            GATK_HF_Command = "{} 12 {}.bam {}".format(GATK_HF_Path, outfile, sset)
                            print(GATK_HF_Command + "\n")
                            os.system(GATK_HF_Command)

                            # GATK_VQSR
                            GATK_VQSR_Command = "{} 12 {}.bam {}".format(GATK_VQSR_Path, outfile, sset)
                            print(GATK_VQSR_Command + "\n")
                            os.system(GATK_VQSR_Command)

                            # Hap.py / rep.py
                            HAPPY_Command = "{}{}.sh".format(HAPPY_prefix, sset)
                            print(HAPPY_Command + "\n")
                            os.system(HAPPY_Command)
                            print("\n")







