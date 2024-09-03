#Mapping: 
bwa mem -t 20 -R '@RG\tID:RYP1ID\tSM:RYP1' index fq1 fq2 | samtools view -@ 20 -bh - | samtools sort -@ 20 -T test -o bam.file - 
#Calling for individual sample:
gatk HaplotypeCaller -R ref.fa –emit-ref-confidence GVCF -I bam.file -O GVCF.file
#Megre GVCF:
gatk GenotypeGVCFs -R ref.fa -V GVCF.vcf ... -O vcf.file
#Annotation:
table_annovar.pl vcf.file Database.dir/ -buildver RYP1  -out prefix -remove -protocol refGene -operation g -nastring . -vcfinput
