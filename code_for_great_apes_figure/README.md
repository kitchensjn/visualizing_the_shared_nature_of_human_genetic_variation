# Steps for recreating the Great Ape Euler diagram

python merge_ancestral_alleles.py

bgzip data/great_apes/ancestral_alleles/genome_chr_pos_aa_GAGP18.txt
tabix -s1 -b2 -e2 data/great_apes/ancestral_alleles/genome_chr_pos_aa_GAGP18.txt.gz

bcftools merge -0 data/great_apes/vcfs/Homo.vcf.gz data/great_apes/vcfs/Pan_paniscus.vcf.gz data/great_apes/vcfs/Pan_troglodytes.vcf.gz -Ou | bcftools view -m2 -M2 -v snps -i 'F_MISSING=0'  -Ou | bcftools annotate -a data/great_apes/ancestral_alleles/genome_chr_pos_aa_GAGP18.txt.gz -h data/great_apes/ancestral_alleles/GAGP_18_genome_header.txt -c CHROM,POS,GAGP_AA_18 -Ou | bcftools view -e "INFO/GAGP_AA_18 = '.'" -Oz -o data/great_apes/vcfs/hbc_genome_biallelic_clean_cleanGAGP_AA.vcf.gz

python hbc_genome_vcfto_freq.py