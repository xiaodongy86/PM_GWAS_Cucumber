
#filter raw seq data, fastp，v0.20.1，> 30% of bases with a sequencing quality score < 20
fastp -i CRR121626_R1.fq.gz -I CRR121626_R2.fq.gz -o CRR121626_f1.fq.gz -O CRR121626_r2.fq.gz--qualified_quality_phred 20 --unqualified_percent_limit 30

#build index, bwa：0.7.17-r1188
bwa index Gy14_genome_v2.fa

#make alignment, generate .sam file
bwa mem Gy14_genome_v2.fa CRR121626_f1.fq.gz CRR121626_r2.fq.gz > CRR121626.sam

# .sam to .bam file
samtools import Gy14_genome_v2.fa CRR121626.sam CRR121626.bam

#Sorting and Indexing
samtools sort CRR121626.bam CRR121626.sorted

# gatk V4.1.4.1
gatk --java-options "-Xmx50g" MarkDuplicates --TMP_DIR tmp --INPUT CRR121626.sorted.bam -OUTPUT CRR121626.markdup.bam \
	--METRICS_FILE CRR121626.markdup.mat --CREATE_INDEX TRUE --REMOVE_DUPLICATES TRUE 


#HaplotypeCaller
gatk HaplotypeCaller -R Gy14_genome_v2.fa -I CRR121626.markdup.bam \
	CRR121626.gatk.raw.g.vcf.gz --emit-ref-confidence GVCF -stand-call-conf 30.0 --min-base-quality-score 10 \
	--dont-use-soft-clipped-bases true --tmp-dir tmp

#CombineGVCFs
gatk CombineGVCFs -R Gy14_genome_v2.fa -V CRR121626.gatk.raw.g.vcf.gz -V CRR121627.gatk.raw.g.vcf.gz \
	-O samples.g.vcf.gz --create-output-variant-index --tmp-dir tmp

#GenotypeGVCFs
gatk GenotypeGVCFs -R Gy14_genome_v2.fa -V samples.g.vcf.gz samples.raw.vcf.gz --tmp-dir tmp

#SelectVariants
gatk --java-options "-Xmx4g" SelectVariants -R Gy14_genome_v2.fa -V samples.gatk.raw.vcf.gz 
	--select-type-to-include SNP --exclude-non-variants -O samples.gatk.con.snp.vcf.gz
gatk --java-options "-Xmx4g" SelectVariants -R Gy14_genome_v2.fa -V samples.gatk.raw.vcf.gz 
	--select-type-to-include INDEL --exclude-non-variants -O samples.gatk.con.indel.vcf.gz

#VQSR
gatk --java-options "-Xmx16g -Xms16g" VariantRecalibrator \
	-R Gy14_genome_v2.fa \
	-V samples.gatk.con.snp.vcf.gz \
	-O samples.gatk.con.snp.recal \
	--tranches-file samples.gatk.con.snp.tranches \
	--rscript-file samples.gatk.con.snp.plots.R \
	-mode SNP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	--resource:hapmap,known=false,training=true,truth=true,prior=10.0 samples.gatk.con.snp.vcf.gz \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 97.0 -tranche 95.0 -tranche 93.0 -tranche 90.0 \
	--max-gaussians 6 --minimum-bad-variants 1000 --bad-lod-score-cutoff -5 \
	--tmp-dir tmp 

gatk ApplyVQSR -mode SNP \
	--truth-sensitivity-filter-level 99 \
	-R Gy14_genome_v2.fa \
	-V samples.gatk.con.snp.vcf.gz \
	--recal-file samples.gatk.con.snp.recal \
	--tranches-file samples.gatk.con.snp.tranches \
	-O samples.gatk.con.snp.vqsr.vcf.gz \
	--tmp-dir tmp

gatk --java-options "-Xmx4g -Xms4g" VariantRecalibrator \
	-R Gy14_genome_v2.fa \
	-V samples.gatk.con.indel.vcf.gz \
	-O samples.gatk.con.indel.recal \
	--tranches-file samples.gatk.con.indel.tranches \
	--rscript-file samples.gatk.con.indel.plots.R \
	-mode INDEL -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	--resource:hapmap,known=false,training=true,truth=true,prior=10.0 samples.gatk.con.indel.vcf.gz \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 97.0 -tranche 95.0 -tranche 93.0 -tranche 90.0 \
	--max-gaussians 6 --minimum-bad-variants 1000 --bad-lod-score-cutoff -5 \
	--tmp-dir tmp

gatk ApplyVQSR -mode INDEL \
	--truth-sensitivity-filter-level 99 \
	-R Gy14_genome_v2.fa \
	-V samples.gatk.con.indel.vcf.gz \
	--recal-file samples.gatk.con.indel.recal \
	--tranches-file samples.gatk.con.indel.tranches \
	-O samples.gatk.con.indel.vqsr.vcf.gz \
	--tmp-dir tmp


gatk VariantFiltration -R Gy14_genome_v2.fa \
		-V samples.gatk.con.snp.vqsr.vcf \
		-O samples.gatk.con.snp.vqsr.filter.vcf \
		--filter-expression "QD < 2.0" --filter-name "QD" \
		--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum" \
		--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum" \
		--filter-expression "FS > 60.0" --filter-name "FS" \
		--filter-expression "MQ < 40.0" --filter-name "MQ" \
		--filter-expression "DP < 4.0" --filter-name "DP" \
		--filter-expression "QUAL < 30.0" --filter-name "QUAL" \
		--tmp-dir tmp
gatk VariantFiltration -R Gy14_genome_v2.fa \
		-V samples.gatk.con.indel.vqsr.vcf \
		-O samples.gatk.con.indel.vqsr.filter.vcf \
		--filter-expression "QD < 2.0" --filter-name "QD" \
		--filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum" \
		--filter-expression "FS > 200.0" --filter-name "FS" \
		--filter-expression "DP < 4.0" --filter-name "DP" \
		--filter-expression "QUAL < 30.0" --filter-name "QUAL" \
		--tmp-dir tmp

perl VarFiltrationVcfs.pl -i samples.gatk.con.snp.vqsr.filter.vcf -o samples.gatk.con.snp.vqsr.filter.result.vcf
perl VarFiltrationVcfs.pl -i samples.gatk.con.indel.vqsr.filter.vcf -o samples.gatk.con.indel.vqsr.filter.result.vcf

#Add id
bcftools annotate --set-id +'%CHROM\_%POS' samples.gatk.con.snp.vqsr.filter.result.vcf -Oz -o test.snp.id.vcf
bcftools annotate --set-id +'%CHROM\_%POS' samples.gatk.con.indel.vqsr.filter.result.vcf -Oz -o test.indel.id.vcf

#filter snp and indel
plink --vcf test.snp.id.vcf --maf 0.05 --geno 0.2 --recode vcf-iid --out test.snp.0.05-0.2.filter --allow-extra-chr
plink --vcf test.indel.id.vcf --maf 0.05 --geno 0.2 --recode vcf-iid --out test.indel.0.05-0.2.filter --allow-extra-chr

#combine 
bcftools sort test.indel.0.05-0.2.filter.vcf -O z -o INDEL_filtered_sorted.vcf.gz
bcftools sort test.snp.0.05-0.2.filter.vcf -O z -o SNP_filtered_sorted.vcf.gz

#index htslib(V:1.9)
tabix -p vcf SNP_filtered_sorted.vcf.gz
tabix -p vcf INDEL_filtered_sorted.vcf.gz

# Merge SNP and INDEL
bcftools concat SNP_filtered_sorted.vcf.gz INDEL_filtered_sorted.vcf.gz  -a -O z -o ALL_filtered_sorted.vcf.gz

#GCTA kinship
plink --vcf ALL_filtered_sorted.vcf --make-bed --out All

gcta64 --bfile All --make-grm-gz --make-grm-alg 1 --out kinship  --autosome-num 19

#plink PCA
vcftools --vcf ALL_filtered_sorted.vcf --plink-tped --out ALL_filtered_pc
plink --noweb --tfile ALL_filtered_pc --make-bed --out ALL_filtered_pc
plink --bfile ALL_filtered_pc --pca 5 --out pca

#Population structure
#filter
plink --vcf ALL_filtered_sorted.vcf --indep-pairwise 100 50 0.2 --out ALL_filtered_sorted --allow-extra-chr
#Output
plink --vcf ALL_filtered_sorted.vcf --make-bed --extract ALL_filtered_sorted.prune.in --out ALL_filtered_sorted.prune.in --allow-extra-chr
#vcf format
plink --bfile ALL_filtered_sorted.prune.in --recode vcf-iid --out ALL_filtered_sorted.prune.in --allow-extra-chr 


#calculate Q
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13;do ./admixture --cv ALL_filtered_ad.bed $K|tee log${K}.out;done

#calculate FST
vcftools --vcf combine.vcf --weir-fst-pop pop_1.txt --weir-fst-pop pop_2.txt --out pop1_vs_pop2
vcftools --vcf combine.vcf --weir-fst-pop pop_2.txt --weir-fst-pop pop_3.txt --out pop2_vs_pop3
vcftools --vcf combine.vcf --weir-fst-pop pop_4.txt --weir-fst-pop pop_1.txt --out pop4_vs_pop1

#Calculate genetic distance
plink --vcf test.filter.prune.in.vcf --distance 1-ibs --out test_filter_prune_in --allow-extra-chr 

#LD decay
PopLDdecay -InVCF ALL_filtered_sorted.vcf -OutStat LDdecay -MaxDist 200 -OutType 2
perl PopLDdecay/bin/Plot_OnePop.pl -inFile LDdecay.stat.gz -output Fig

# -- GWAS TASSEL\ 5/
#vcf to hapmap
perl run_pipeline.pl -Xms512m -Xmx12g -vcf ALL_filtered_sorted.vcf -export work -exportType Hapmap

#sort SNP and InDel sites
perl run_pipeline.pl -Xms512m -Xmx12g  -SortGenotypeFilePlugin -inputFile test.hmp.txt -outputFile test_sort -fileType Hapmap 

#GLM
perl run_pipeline.pl -Xms512m -Xmx12g -fork1 -h test_sort.hmp.txt -fork2 -r trait02.txt -fork3 -q Q.txt \
	-excludeLastTrait -combine4 -input1 -input2 -input3 -intersect -glm \
	-export test_glm02 -runfork1 -runfork2 -runfork3

#MLM
perl run_pipeline.pl -Xms512m -Xmx12g -fork1 -h test_sort.hmp.txt -fork2 -r trait02.txt -fork3 -q Q.txt \
	-excludeLastTrait -fork4 -k K.txt -combine5 -input1 -input2 -input3 -intersect -combine6 -input5 -input4 \
	-mlm -mlmVarCompEst P3D -mlmCompressionLevel None \
	-export test_mlm02 -runfork1 -runfork2 -runfork3 -runfork4

