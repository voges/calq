import os.path

# Parameters
datasets = [["ERP001775", "ERR174324.aln_bowtie2.sorted.dupmark.rg.realn.recal"],
            ["NA12878_Garvan_replicate_J", "NA12878_V2.5_Robot_2.aln_bowtie2.sorted.dupmark.rg.realn.recal"],
            ["NA12878-SRX517292", "SRR1238539.aln_bowtie2.sorted.dupmark.rg.realn.recal"]]
subsets = ["3", "11", "20"]
filtersize = ["10", "17", "25"]
filtertype = ["Gauss"]
quantizerType = ["Uniform"]
quantSteps = [["2", "8"]]
squashed = [""]


codecs = [[".sam", ".calq-0b74e3d", ""], [".bam", ".crumble-0.5", "-1"],[".bam", ".crumble-0.5", "-9"], ["", ".dsrc_2.0-il8b", ""], [".fastq", ".quartz-0.2.2", ".fastq.clipped_qual.aln_bowtie2.sorted.dupmark.rg"], [".sam", ".qvz2-d5383c6", "-t1"], [".sam", ".qvz2-d5383c6", "-t2"], ["", "", ""], [".sam", ".qvz2-d5383c6", "-t4"], [".sam", ".qvz2-d5383c6", "-t8"]]

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
        for codec in codecs:
            outfile = "{}{}/{}.{}{}{}{}".format(folder, codec[1],  dset[1], sset, codec[0], codec[1], codec[2])

            if not os.path.isfile(outfile + ".Platypus.log"):
                # Platypus
                PlatypusCommand = "{} 12 {}.bam {}".format(platypusPath, outfile, sset)
                print(PlatypusCommand + "\n", flush=True)
                os.system(PlatypusCommand)
                os.system("mv .bam.GATK_HF.log {}.Platypus.log".format(outfile))
            else:
                print("{}.Platypus.log already existing. Skipping platypus.".format(outfile) + "\n", flush=True)

            if not os.path.isfile(outfile + ".GATK_HF.log"):
                # GATK_HF
                GATK_HF_Command = "{} 12 {}.bam {}".format(GATK_HF_Path, outfile, sset)
                print(GATK_HF_Command + "\n", flush=True)
                os.system(GATK_HF_Command)
                os.system("mv .bam.GATK_HF.log {}.GATK_HF.log".format(outfile))
            else:
                print("{}.GATK_HF.log already existing. Skipping GATK_HF.".format(outfile) + "\n", flush=True)

            if not os.path.isfile(outfile + ".GATK_VQSR.log"):
                # GATK_VQSR
                GATK_VQSR_Command = "{} 12 {}.bam {}".format(GATK_VQSR_Path, outfile, sset)
                print(GATK_VQSR_Command + "\n", flush=True)
                os.system(GATK_VQSR_Command)
                os.system("mv .bam.GATK_VQSR.log {}.GATK_VQSR.log".format(outfile))
            else:
                print("{}.GATK_VQSR.log already existing. Skipping GATK_VQSR.".format(outfile) + "\n", flush=True)







