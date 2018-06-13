import os.path
import csv

# Parameters
datasets = [["ERP001775", "ERR174324.aln_bowtie2.sorted.dupmark.rg.realn.recal"],
            ["NA12878_Garvan_replicate_J", "NA12878_V2.5_Robot_2.aln_bowtie2.sorted.dupmark.rg.realn.recal"],
            ["NA12878-SRX517292", "SRR1238539.aln_bowtie2.sorted.dupmark.rg.realn.recal"]]
subsets = ["20"]
filtersize = ["17"]
filtertype = ["Gauss"]
quantizerType = ["Uniform"]
quantSteps = [["2", "8"]]
squashed = [""]


basedir = "/data/voges/muenteferi"
install_path = "/project/dna/install"
samtools = install_path + "/samtools-1.3/bin/samtools"
calqPath = "/home/muenteferi/Dokumente/calqBuild/calq"
scramble = "/project/dna/bin/scramble"
cramsize = "/project/dna/bin/cram_size"
referencePath = "/data/voges/muenteferi/GATK_bundle-2.8-b37/human_g1k_v37.fasta"
replacePath = "/home/muenteferi/Dokumente/calq/src/ngstools/replace_qual_sam.py"
platypusPath = "/home/muenteferi/Dokumente/calq/src/variant_calling_pipelines/Platypus.sh"
GATK_HF_Path = "/home/muenteferi/Dokumente/calq/src/variant_calling_pipelines/GATK_HF.sh"
GATK_VQSR_Path = "/home/muenteferi/Dokumente/calq/src/variant_calling_pipelines/GATK_VQSR.sh"
HAPPY_prefix = "/home/muenteferi/Dokumente/calq/src/variant_calling_pipelines/do_happy+reppy_chr"
vcfList = [".platypus.snps.vcf", ".GATK.snps.hard_filtered.vcf", ".GATK.snps.filtered900.vcf",
           ".GATK.snps.filtered990.vcf", ".GATK.snps.filtered999.vcf", ".GATK.snps.filtered1000.vcf",
           ".platypus.indels.vcf", ".GATK.indels.hard_filtered.vcf", ".GATK.indels.filtered900.vcf",
           ".GATK.indels.filtered990.vcf", ".GATK.indels.filtered999.vcf", ".GATK.indels.filtered1000.vcf"]
outCSV = "results.csv"

if os.path.isfile(outCSV):
    print("Outfile already exists. Exiting.",flush=True)
    exit(-1)

f = open(outCSV, 'w')
f.write("Data Chromosome Codec Caller Truth TruePositive FalseNegative FalsePositive Unknown precision recall f-score Filesize\n");

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
                            filesize = os.path.getsize(outfile)

                            # Hap.py / rep.py
                            for vcf in vcfList:
                                if not os.path.isfile("{}.bam".format(outfile)):
                                    print("File '{}.bam' does not exist. Skipping!\n".format(outfile), flush=True)
                                    continue

                                happyCSV = "{}.bam{}.happy.summary.csv".format(outfile, vcf)
                                HAPPY_Command = "{}{}.sh {}.bam{}".format(HAPPY_prefix, sset, outfile, vcf)

                                if os.path.isfile(happyCSV):
                                    print("File {} already existing. Skipping Happy, just collecting info!\n".format(happyCSV), flush=True)
                                else:
                                    print(HAPPY_Command + "\n", flush=True)
                                    os.system(HAPPY_Command)

                                precision = 0.0
                                recall = 0.0
                                fScore = 0.0
                                calls = 0
                                truePositive = 0
                                falsePositive = 0
                                unknownCalls = 0
                                falseNegative = 0
                                rowCtr = 0
                                colCtr = 0

                                if not os.path.isfile(happyCSV):
                                    print("FAIL!\n", flush=True)
                                    f.write("{} {} calq-haplo{}{}{} {} {} {} {} {} {} {} {} {} {}\n".format(dset[0], sset, fsize, qtype, squash, vcf, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, filesize))
                                    continue
                                    

                                with open(happyCSV, newline='') as File:
                                    reader = csv.reader(File)
                                    for row in reader:
                                        for col in row:
                                            if ".snps." in vcf:
                                                if rowCtr == 4 and colCtr == 9:
                                                    recall = float(col)
                                                if rowCtr == 4 and colCtr == 10:
                                                    precision = float(col)
                                                if rowCtr == 4 and colCtr == 2:
                                                    calls = int(col)
                                                if rowCtr == 4 and colCtr == 3:
                                                    truePositive = int(col)
                                                if rowCtr == 4 and colCtr == 4:
                                                    falseNegative = int(col)
                                                if rowCtr == 4 and colCtr == 6:
                                                    falsePositive = int(col)
                                                if rowCtr == 4 and colCtr == 7:
                                                    unknownCalls = int(col)

                                            else:
                                                if rowCtr == 2 and colCtr == 9:
                                                    recall = float(col)
                                                if rowCtr == 2 and colCtr == 10:
                                                    precision = float(col)
                                                if rowCtr == 2 and colCtr == 2:
                                                    calls = int(col)
                                                if rowCtr == 2 and colCtr == 3:
                                                    truePositive = int(col)
                                                if rowCtr == 2 and colCtr == 4:
                                                    falseNegative = int(col)
                                                if rowCtr == 2 and colCtr == 6:
                                                    falsePositive = int(col)
                                                if rowCtr == 2 and colCtr == 7:
                                                    unknownCalls = int(col)


                                            colCtr = colCtr + 1
                                        rowCtr = rowCtr + 1
                                        colCtr = 0
                                    fScore = 0.0 if recall == 0.0 and precision == 0.0 else 2.0 * precision * recall / (precision + recall)
                                    print("P: {} R: {} F: {};\n".format(precision, recall, fScore), flush=True)
                                    f.write("{} {} calq-haplo{}{}{} {} {} {} {} {} {} {} {} {} {}\n".format(dset[0], sset, fsize, qtype, squash, vcf, calls, truePositive, falseNegative, falsePositive, unknownCalls, precision, recall, fScore, filesize))
                                print("\n", flush=True)
f.close()








