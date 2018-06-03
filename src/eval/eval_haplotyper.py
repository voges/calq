import os.path

# Parameters
datasets = [["ERP001775", "ERR174324.aln_bowtie2.sorted.dupmark.rg.realn.recal"],
            ["NA12878_Garvan_replicate_J", "NA12878_V2.5_Robot_2.aln_bowtie2.sorted.dupmark.rg.realn.recal"],
            ["NA12878-SRX517292", "SRR1238539.aln_bowtie2.sorted.dupmark.rg.realn.recal"]]
subsets = ["3", "11", "20"]
filtersize = ["10"]
filtertype = ["Gauss"]
quantizerType = ["Uniform"]
quantSteps = [["2", "8"]]
squashed = [""]

# Paths
basedir = "/data/voges/muenteferi"
install_path = "/project/dna/install"
samtools = install_path + "/samtools-1.3/bin/samtools"
calqPath = "/home/muenteferi/Dokumente/calqBuild/calq"
referencePath = "/data/voges/muenteferi/GATK_bundle-2.8-b37/human_g1k_v37.fasta"
replacePath = "/home/muenteferi/Dokumente/calq/src/ngstools/replace_qual_sam.py"
platypusPath = "/home/muenteferi/Dokumente/calq/src/variant_calling_pipelines/Platypus.sh"
GATK_HF_Path = "/home/muenteferi/Dokumente/calq/src/variant_calling_pipelines/GATK_HF.sh"
GATK_VQSR_Path = "/home/muenteferi/Dokumente/calq/src/variant_calling_pipelines/GATK_VQSR.sh"
HAPPY_prefix = "/home/muenteferi/Dokumente/calq/src/variant_calling_pipelines/do_happy+reppy_chr"

vcfList = [".platypus.snps.vcf", ".GATK.snps.hard_filtered.vcf", ".GATK.snps.filtered900.vcf",
           ".GATK.snps.filtered990.vcf", ".GATK.snps.filtered999.vcf", ".GATK.snps.filtered1000.vcf"]

for dset in datasets:
    for sset in subsets:
        folder = "{}/{}/{}.{}".format(basedir, dset[0], dset[1], sset)
        filename = "{}.{}".format(dset[1], sset)
        filepath = "{}/{}".format(folder, filename)
        if not os.path.isfile(filepath + ".bam"):
            print("NOT existing: " + filepath + ".bam", flush=True)
            exit(-1)

        if not os.path.isfile(filepath + ".sam"):
            print("Converting to sam: " + file + ".bam" , flush=True)
            command = "{} view -h -o {}.sam {}.bam".format(samtools, filepath, filepath)
            os.system(command)

        for fsize in filtersize:
            for ftype in filtertype:
                for qtype in quantizerType:
                    for qsteps in quantSteps:
                        for squash in squashed:
                            # Create folders and SAM input
                            outfolder = "{}.calq-haplo.filter{}{}.quant{}_{}{}{}".format(folder,fsize,ftype,qsteps[0], qsteps[1], qtype, squash)
                            outfile = "{}.calq-haplo.filter{}{}.quant{}_{}{}{}.cq".format(outfolder + "/" + filename, fsize, ftype, qsteps[0], qsteps[1], qtype, squash)
                            if not os.path.isdir(outfolder):
                                print("Creating dir: " + outfolder, flush=True)
                                os.system("mkdir " + outfolder)

                            # Compress
                            calqCommand = "{} -q Illumina-1.8+ -p 2 -b 10000 {}.sam  -o {} -r {} --quantizerType {} " \
                                          "--filterType {} --quantizationMin {} --quantizationMax {} --filterSize {} {}".\
                                format(calqPath, filepath, outfile, referencePath, qtype, ftype, qsteps[0],
                                       qsteps[1], fsize, squash)
                            print(calqCommand + "\n", flush=True)
                            os.system(calqCommand)

                            # Decompress
                            calqDeCommand = "{} -d -s {}.sam {} -o {}.qual".format(calqPath, filepath, outfile, outfile)
                            print(calqDeCommand + "\n", flush=True)
                            os.system(calqDeCommand)

                            # Insert quality values
                            replaceCommand = "python {} {}.sam {}.qual 1> {}.sam".format(replacePath, filepath, outfile, outfile)
                            print (replaceCommand + "\n", flush=True)
                            os.system(replaceCommand)

                            # Create bam
                            convertCommand = "{} view -bh {}.sam > {}.bam".format(samtools, outfile, outfile)
                            print(convertCommand + "\n", flush=True)
                            os.system(convertCommand)

                            # Index bam
                            indexCommand = "{} index -b {}.bam {}.bai".format(samtools, outfile, outfile)
                            print(indexCommand + "\n", flush=True)
                            os.system(indexCommand)

                            # Platypus
                            os.system("rm .bam.GATK_HF.log")
                            PlatypusCommand = "{} 12 {}.bam {}".format(platypusPath, outfile, sset)
                            print(PlatypusCommand + "\n", flush=True)
                            os.system(PlatypusCommand)
                            os.system("mv .bam.GATK_HF.log {}.Playtypus.log".format(outfile))


                            # GATK_HF
                            os.system("rm .bam.GATK_HF.log")
                            GATK_HF_Command = "{} 12 {}.bam {}".format(GATK_HF_Path, outfile, sset)
                            print(GATK_HF_Command + "\n", flush=True)
                            os.system(GATK_HF_Command)
                            os.system("mv .bam.GATK_HF.log {}.GATK_HF.log".format(outfile))

                            # GATK_VQSR
                            os.system("rm .bam.GATK_VQSR.log")
                            GATK_VQSR_Command = "{} 12 {}.bam {}".format(GATK_VQSR_Path, outfile, sset)
                            print(GATK_VQSR_Command + "\n", flush=True)
                            os.system(GATK_VQSR_Command)
                            os.system("mv .bam.GATK_VQSR.log {}.GATK_VQSR.log".format(outfile))







