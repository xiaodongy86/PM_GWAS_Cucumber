minimap2 -ax asm5 -t 90 -x intractg -t 80 ref.fa contigs.fa | smatools view -bh - > bam.file
python assemblatron.py --sv --bam bam.file > vcf1.file
java -jar npInv1.26.jar --input bam.file --output vcf2.file --min 50 --max 10000000
pbmm2 align ref.fa input.fa bam.file --sort --preset CCS  --sample prefix
pbsv discover bam.file svsig.gz
pbsv call --ccs ref.fa svsig.gz vcf3.file
SURVIVOR merge SV.file.list 100 1 1 1 0 30 merged.vcf
