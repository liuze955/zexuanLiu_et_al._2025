Both read alignment and variant calling were performed by GTX-One variants caller. The code below shows an example.
```shell
# fastq align
gtx wgs -g -o ind1.g.vcf.gz -R '@RG\tID:ind1\tSM:ind1' GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna ind1.fastq1 ind1.fastq2 -t 48
```
Then, perform a joint calling
```shell
# Import GVCFs into Genomics DB database
gtx gi \
  -v data/gvcfs/ind1.g.vcf.gz \
  -v data/gvcfs/ind2.g.vcf.gz \
  ...
  --genomicsdb-workspace-path my_database
  -r GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna

# joint calling
gtx genotype_gvcfs \
   -r GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna \
   -v gendb://my_database \
   -o output.vcf.gz
```
Filter vcf files
```shell
# gatk hardfilter
gatk SelectVariants -R GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna -V output.vcf.gz -select-type SNP -O snp.vcf.gz

gatk VariantFiltration -R GCF_016772045.1_ARS-UI_Ramb_v2.0_genomic.fna -V snp.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30"  -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -O snps_hardFiltered.vcf.gz

vcftools --gzvcf snps_hardFiltered.vcf.gz --remove-filtered-all --recode --out snps_filter_pass

# Retain biallelic sites
vcftools --vcf snps_filter_pass.recode.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out snp_biallelic

# Change chromosome names to numbers
bcftools annotate --rename-chrs old.new.name snp_biallelic.recode.vcf  -Oz -o sheep.newchr.vcf.gz --threads 10

# Add ID for plink filtering
bcftools annotate --set-id +'%CHROM\_%POS' sheep.newchr.vcf.gz  -Ov -o id.vcf  --threads 10

# SNP missing rate filtering
plink --vcf id.vcf --geno 0.1 --keep-allele-order --allow-extra-chr --sheep  --recode  vcf-iid --out genofilter

# individual missing rate filtering
plink --vcf genofilter.vcf --mind 0.2 --keep-allele-order --allow-extra-chr --sheep  --recode  vcf-iid --out mindfilter

# phase using beagle
java -Xmx100g -jar beagle.08Feb22.fa4.jar gt=genofilter.vcf out=beagle nthreads=28

# maf filtering
vcftools --gzvcf beagle.vcf.gz --maf 0.01 --recode --recode-INFO-all --out maf

# LD filtering
plink --vcf maf.recode.vcf --allow-extra-chr --keep-allele-order --indep-pairwise 50 10 0.2 --sheep  --out ld
vcftools --vcf maf.recode.vcf --snps ld.prune.in --recode --recode-INFO-all --out article
```
In order to avoid the influence of sample size on population structure, we sampled Hu sheep, and the retained samples are shown in equal.sample.txt
```R
result <- info[info$sort == "Hu",] %>%
  group_by(pop) %>%
  sample_n(size = 2)
```

```shell
awk '{print $1,$1}' equal.sample.txt >　ind.keep
plink --bfile article --keep ind.keep --sheep --keep-allele-order --make-bed --recode vcf-iid --out equal
```
