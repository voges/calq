import os.path
import csv

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


codecs = [[".sam", ".calq-0b74e3d", ""], [".bam", ".crumble-0.5", "-1"],[".bam", ".crumble-0.5", "-9"], ["", ".dsrc-2.0-il8b", ""], ["", ".quartz-0.2.2", ".clipped_qual.aln_bowtie2.sorted.dupmark.rg"], [".sam", ".qvz2-d5383c6", "-t1"], [".sam", ".qvz2-d5383c6", "-t2"], ["", "", ""], [".sam", ".qvz2-d5383c6", "-t4"], [".sam", ".qvz2-d5383c6", "-t8"]]

# Paths
basedir = "/data/voges/muenteferi"
install_path = "/project/dna/install"
samtools = install_path + "/samtools-1.3/bin/samtools"
scramble = "/project/dna/bin/scramble"
cramsize = "/project/dna/bin/cram_size"
calqPath = "/home/muenteferi/Dokumente/calqBuild/calq"
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

outCSV = "resultsCodec.csv"
f = open(outCSV, 'w')
f.write("Data Chromosome Codec Caller Truth TruePositive FalseNegative FalsePositive Unknown precision recall f-score Filesize\n");

for dset in datasets:
    for sset in subsets:
        folder = "{}/{}/{}.{}".format(basedir, dset[0], dset[1], sset)
        filename = "{}.{}".format(dset[1], sset)
        filepath = "{}/{}".format(folder, filename)
        for codec in codecs:
            outfile = "{}{}/{}.{}{}{}{}".format(folder, codec[1],  dset[1], sset, codec[0], codec[1], codec[2])
            cname = codec[1]+codec[2]
            if cname == "":
                cname = ".none"
            if cname == ".quartz-0.2.2.clipped_qual.aln_bowtie2.sorted.dupmark.rg":
                cname = ".quartz-0.2.2"

            filesize = "0"

            if not os.path.isfile(outfile + ".bam.cram"):
                CRAM_Command = "{} -r {} {}.bam {}.bam.cram".format(scramble, referencePath, outfile, outfile)
                os.system(CRAM_Command)
            else:
                print("File '{}.bam.cram' exists. Skipping!\n".format(outfile), flush=True)

            Size_Command = "{} {}.bam.cram | grep QS | awk \'{print $6}\'".format(cramsize, outfile)
            filesize = os.popen(Size_Command).read()

            # Hap.py / rep.py
            for vcf in vcfList:
                if not os.path.isfile("{}.bam".format(outfile)):
                    print("File '{}.bam' does not exist. Skipping!\n".format(outfile), flush=True)
                    continue

                happyCSV = "{}.bam{}.happy.summary.csv".format(outfile, vcf)
                HAPPY_Command = "{}{}.sh {}.bam{}".format(HAPPY_prefix, sset, outfile, vcf)

                if os.path.isfile(happyCSV):
                    print("File {} already existing. Skipping Happy, just collecting info!\n".format(happyCSV),
                          flush=True)
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
                    f.write("{} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(dset[0], sset, cname, vcf, -1.0, -1.0,
                                                                                            -1.0, -1.0, -1.0, -1.0,
                                                                                            -1.0, -1.0, -1.0))
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
                    fScore = 0.0 if recall == 0.0 and precision == 0.0 else 2.0 * precision * recall / (
                                precision + recall)
                    print("P: {} R: {} F: {};\n".format(precision, recall, fScore), flush=True)
                    f.write("{} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(dset[0], sset, cname, vcf, calls,
                                                                                            truePositive, falseNegative,
                                                                                            falsePositive, unknownCalls,
                                                                                            precision, recall, fScore,
                                                                                            filesize))
                print("\n", flush=True)
f.close()








