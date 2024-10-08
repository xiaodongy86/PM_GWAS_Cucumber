#Mapping
#Calling for individual sample:
python configManta.py –bam bam.file --runDir target.dir --referenceFasta ref.fa
python runWorkflow.py
#Merging:
SURVIVOR merge SV.file.list 100 1 1 1 0 30 merged.vcf
#Genotyping:
idxdepth -b bam.file -o out.dir -r ref.fa --threads 40
multigrmpy.py --threads 20 -i merged.vcf-m config.file -r ref.fa -o out.dir
