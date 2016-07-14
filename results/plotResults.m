% the files should be in CSV format without header line(s)
%fileName = '/project/dna/git/calq/test_files/DH10B+MiSeq_Ecoli_DH10B_110721_PF.sam.200000Lines.sam.csv';
%fileName = '/project/dna/git/calq/test_files/test.sam.csv';
%fileName = '/data/gidb/NA12878/tmp/NA12878.pacbio.bwa-sw.20140202.sam.csv';
fileName = '/project/dna/git/calq/build_linux_gcc/tmp';

M = csvread(fileName);

figure(1);
%plot(M(:,1)); % locus

plot(M(:,2)); % sequencing depth
hold on;

plot(M(:,3)); % confidence
%plot(M(:,4)); % quantizer index

