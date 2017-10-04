#!/bin/bash

happy_reppy_sh="/home/voges/git/calq/src/variant_calling_pipelines/happy+reppy.sh"
num_threads=4
variants_vcf=$1
log_txt="$variants_vcf.happy+reppy.log"
ref_fasta="/data/gidb/tmp/MPEG/GATK_bundle-2.8-b37/human_g1k_v37.fasta"
golden_vcf_gz="/data/gidb/tmp/MPEG/GIAB-NA12878_HB001-NISTv3.2.2.20/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.20.vcf.gz"
golden_bed="/data/gidb/tmp/MPEG/GIAB-NA12878_HB001-NISTv3.2.2.20/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.20.bed"

bash $happy_reppy_sh $num_threads $log_txt $variants_vcf $ref_fasta $golden_vcf_gz $golden_bed

