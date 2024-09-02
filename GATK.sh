genome=Gy14_genome_v2.fa

GATK -T HaplotypeCaller -R $genome -I bamfile \
			 -o sample.gatk.raw.g.vcf \
			 --emitRefConfidence GVCF -stand_call_conf 30.0 -rf BadCigar --variant_index_type LINEAR \
			 --variant_index_parameter 128000 -nct 5

GATK -T GenotypeGVCFs -R $genome -V sample1.gatk.raw.g.vcf -V sample2.gatk.raw.g.vcf  -o samples.gatk.con.vcf


GATK -T SelectVariants -R $genome -V samples.gatk.con.vcf -selectType SNP -o samples.gatk.con.snp.vcf

GATK -T VariantRecalibrator -R $genome -input samples.gatk.con.snp.vcf \
	-recalFile samples.gatk.con.snp.recal \
	-tranchesFile samples.gatk.con.snp.tranches \
	-rscriptFile samples.gatk.con.snp.plots.R -mode SNP -an DP -an MQRankSum -an ReadPosRankSum \
	-resource:hapmap,known=false,training=true,truth=true,prior=10.0 samples.gatk.con.vcf

GATK -T ApplyRecalibration -mode SNP --ts_filter_level 90 -recalFile samples.gatk.con.snp.recal \
	-tranchesFile samples.gatk.con.snp.tranches \
	-R samples.gatk.con.snp.vcf \
	-o samples.gatk.con.snp.vqsr.vcf

GATK -T VariantFiltration -R $genome \
		-V samples.gatk.con.snp.vqsr.vcf \
		-o samples.gatk.con.snp.vqsr.filter.vcf \
		--filterExpression "QD < 2.0" --filterName "QD" --filterExpression "SOR > 6.0" --filterName "SOR" \
		--filterExpression "QUAL < 5" --filterName "QUAL"
