filterExpression="QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
filterName="GATK_Recommended"

java -jar /home/mikel/Programs/GenomeAnalysisTK.jar -R $pathHumanReference -T VariantFiltration -L $targetRegion -V $f --filterExpression "$filterExpression" --filterName "$filterName" -o $outputFilteredVCF
