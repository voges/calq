%fileName = '/project/dna/git/calq/test_files/DH10B+MiSeq_Ecoli_DH10B_110721_PF.sam.200000Lines.sam.stats';
%fileName = '/project/dna/git/calq/test_files/polyploidyTest.sam.stats';
%fileName = '/project/dna/git/calq/test_files/test.sam.stats';
fileName = '/data/gidb/NA12878/tmp/NA12878.pacbio.bwa-sw.20140202.sam.stats';

fileID = fopen(fileName);
while ~feof(fileID);
    C = textscan(fileID,'%f%s%f%f','Delimiter',',');
end;
fclose(fileID);

figure(1);
plot(C{1,1}); % sequencing depth
hold on;
%plot(-log(C{1,3})); % entropy
plot(C{1,4}); % diff
