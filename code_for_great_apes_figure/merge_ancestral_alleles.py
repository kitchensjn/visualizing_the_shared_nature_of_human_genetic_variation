import glob
import gzip


with open("../data/great_apes/ancestral_alleles/genome_chr_pos_aa_GAGP18.txt", "w") as GAGPlist:
    for chrom in range(1,23):
        with gzip.open(f"../data/great_apes/ancestral_alleles/raw/chr{chrom}.calls.txt.gz", 'rt') as gzippath:
        # read = gzippath.read()
            next(gzippath)
            for line in gzippath:
                list = line.split("\t")
                inner18 = list[19]
                if "?" not in inner18 and "/" not in inner18:
                    GAGPlist.write(list[0] + "\t" + list[1] + '\t' + list[19] + '\n')

with open("data/great_apes/ancestral_alleles/GAGP_18_genome_header.txt", "w") as GAGPheader:
    GAGPheader.write('##INFO=<ID=GAGP_AA_18,Number=1,Type=String,Description="Ancestral allele from GAGP, inner 18, whole genome">')