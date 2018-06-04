import os.path
import csv

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
outCSV = "results.csv"

if os.path.isfile(outCSV):
    print("Outfile already exists. Exiting.",flush=True)
    exit(-1)

f = open(outCSV, 'w')
f.write("Data Chromosome Filtersize Filter Quantizer QuantizerMin QuantizerMax SquashMode Caller IndelP IndelR IndelF SNPP SNPR SNPF\n");

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

                            # Hap.py / rep.py
                            for vcf in vcfList:
                                if not os.path.isfile("{}.bam".format(outfile)):
                                    print("File '{}.bam' does not exist. Skipping!\n".format(outfile), flush=True)
                                    continue

                                happyCSV = "{}.bam{}.happy.summary.csv".format(outfile, vcf)
                                HAPPY_Command = "{}{}.sh {}.bam{}".format(HAPPY_prefix, sset, outfile, vcf)

                                if os.path.isfile(happyCSV):
                                    print("File already existing. Skipping Happy, just collecting info!\n".format(outfile), flush=True)
                                else:
                                    print(HAPPY_Command + "\n", flush=True)
                                    os.system(HAPPY_Command)

                                indelP = 0.0
                                indelR = 0.0
                                SNPP = 0.0
                                SNPR = 0.0
                                rowCtr = 0.0
                                colCtr = 0.0
                                indelF = 0.0
                                SNPF = 0.0

                                with open(happyCSV, newline='') as File:
                                    reader = csv.reader(File)
                                    for row in reader:
                                        for col in row:
                                            if rowCtr == 4 and colCtr == 9:
                                                SNPR = float(col)
                                            if rowCtr == 4 and colCtr == 10:
                                                SNPP = float(col)
                                            if rowCtr == 2 and colCtr == 9:
                                                indelR = float(col)
                                            if rowCtr == 2 and colCtr == 10:
                                                indelP = float(col)
                                            colCtr = colCtr + 1
                                        rowCtr = rowCtr + 1
                                        colCtr = 0
                                    indelF = 0.0 if indelR == 0.0 and indelP == 0.0 else 2.0 * indelP * indelR / (
                                                indelP + indelR)
                                    SNPF = 0.0 if SNPR == 0.0 and SNPP == 0.0 else 2.0 * SNPP * SNPR / (SNPP + SNPR)
                                    print("indel: P: {} R: {} F: {}; SNP: P: {} R: {} F: {}\n".format(indelP, indelR,
                                                                                                    indelF, SNPP, SNPR,
                                                                                                    SNPF), flush=True)
                                    f.write("{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(dset[0], sset, fsize, ftype, qtype, qsteps[0], qsteps[1], squash, vcf, indelP, indelR,
                                                                                                    indelF, SNPP, SNPR,
                                                                                                    SNPF))
                                print("\n", flush=True)
f.close()








