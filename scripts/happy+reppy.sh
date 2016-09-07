#!/bin/bash

###############################################################################
#                Script to report variant calling results with                #
#                 happy (https://github.com/Illumina/hap.py)                  #
#                                    and                                      #
#              reppy (https://github.com/ga4gh/benchmarking-tools/)           #
###############################################################################

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 giab|platinum variants.vcf ref.fasta"
    exit -1
fi

compare=$1
variants_VCF=$2
ref_FASTA=$3

### Ground truth
if [ compare == "giab" ]; then
    ground_truth_VCF_GZ="/data/gidb/simulations/GIAB-NA12878_HB001-NISTv3.2.2/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.vcf.gz"
    ground_truth_BED="/data/gidb/simulations/GIAB-NA12878_HB001-NISTv3.2.2/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed"
elif [ compare == "platinum" ]; then
    ground_truth_VCF_GZ="/data/gidb/simulations/Illumina_Platinum_Genomes-hg19-8.0.1-NA12878/NA12878.vcf.gz"
    ground_truth_BED="/data/gidb/simulations/Illumina_Platinum_Genomes-hg19-8.0.1-NA12878/ConfidentRegions.bed"
else
    echo "Usage: $0 giab|platinum variants.vcf ref.fasta"
    exit -1
fi

### Programs
install_path="/project/dna/install"
hap_py="$install_path/hap.py-0.3.1/bin/hap.py"
rep_py="$install_path/benchmarking-tools-c458561/reporting/basic/bin/rep.py"

### File names
base=$(echo $variants_VCF | sed 's/\.[^.]*$//') # strip .vcf
happy_root=$base".happy"
reppy_HTML=$base".reppy.html"

### Run hap.py
python $hap_py --threads $num_threads --verbose $ground_truth_VCF_GZ $variants_VCF -f $ground_truth_BED -o $happy_root -r $ref_FASTA --roc VQLSOD

### Run rep.py
#rm -f "$happy_root".rep.tsv
printf "method\tcomparisonmethod\tfiles\n" >> "$happy_root".rep.tsv
printf "gatk3\t$compare\t" >> "$happy_root".rep.tsv
printf "$happy_root.extended.csv," >> "$happy_root".rep.tsv
for i in "$happy_root".roc.Locations.*; do printf "$i,"; done >> "$happy_root".rep.tsv
sed -i '$ s/.$//' "$happy_root".rep.tsv
python $rep_py -o $reppy_HTML -l "$happy_root".rep.tsv
