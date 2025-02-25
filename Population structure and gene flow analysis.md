# Population structure 

## PCA
We performed all PCA analysis using GCTA
```shell
gcta --bfile equal --make-grm --thread-num 10 --out geno_grm --autosome-num 26
gcta --grm geno_grm --pca 5 --out out_pca
```
## Admixture
```shell
for K in {2..10}
do
admixture --cv equal.bed $K -j5 |tee log$K.out 
done
```

# Genetic diversity and LD calculation
## Het
```shell
vcftools --vcf article.vcf --het --recode --out het
```
## Pi
```shell
vcftools --vcf article.vcf --out pop.pi --site-pi --keep pop
```
## LD
```shell
LDBlockShow -InVCF article.vcf.gz -OutPut unc5c_bmpr1b -Region 6:29597209:30136603 -OutPng -SeleVar 1 -SpeSNPName snp
```

# Phylogenetic and gene flow analysis
## treemix
```shell
# Prepare input file
plink --bfile article --transpose --sheep --keep-allele-order --recode --out article.t
plink --tfile article.t --freq --within sample.pop.cov --sheep --keep-allele-order

gzip -c plink.frq.strat > plink.frq.strat.gz
python2 plink2treemix.py plink.frq.strat.gz treemix.in.gz

# calculate
for i in {1..10}  #number of introgression
do
        cd work_path
        mkdir result_mn$i
        cd result_mn$i
        for x in {1..5}
        do
            treemix -i treemix.in.gz  -o result_mn$i/m$i.$x -root OUT -k 1000 -bootstrap -global -m $i 
        done
done
```

Select the best result to draw using plotting_funcs.R in treemix
```R
# select the best result
library(OptM)
test = optM("result_m10/")
plot_optM(test)

# draw
source("plotting_funcs.R")   
plot_tree("m1.4")
```


## Fbranch
```shell
Dsuite Dtrios article.vcf set.txt -t treemix.same.txt -o dtrois.result
Dsuite Fbranch treemix.same.txt dtrois.result_tree.txt > fbranch.out
```

## qpGraph

```R
library(admixtools)
library(tidyverse)
library(openxlsx)

my_f2_dir = 'D:/Work/R/highRseqSheep24/admix2/f2/'
f2_blocks = f2_from_precomp(my_f2_dir,afprod = TRUE)

# find graph,numadmix = 2 as an example
opt_results = find_graphs(f2_blocks, numadmix = 2, outpop = 'OUT',
                          stop_gen = 100,seed=123)
winner = opt_results %>% slice_min(score, with_ties = FALSE)
winner$score[[1]]
plot_graph(winner$edges[[1]],color = F)
write.table(winner$edges[[1]],"topo.num2.best.txt",quote = F,row.names = F)

# out of sample score
nblocks = dim(f2_blocks)[3]
train = sample(1:nblocks, round(nblocks/2))
res = qpgraph(data = f2_blocks[,,train], graph,
              f2_blocks_test = f2_blocks[,,-train])
res$score
```

# Selection signature
## Fst
```shell
vcftools --vcf article.vcf --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --out pop1_pop2.txt --fst-window-size 50000 --fst-window-step 25000
```
## Generalized mixed linear model
```shell
# step1. phe中湖羊是1，单羔群体是0
data_path='article.hu_single'
Rscript step1_fitNULLGLMM.R \
        --plinkFile=$data_path \
        --phenoFile=phe \
        --phenoCol=y_binary \
        --sampleIDColinphenoFile=IID \
        --traitType=binary \
        --outputPrefix=step1 \
        --nThreads=24 \
        --IsOverwriteVarianceRatioFile=TRUE \
        --includeNonautoMarkersforVarRatio=TRUE

# step2
Rscript step2_SPAtests.R \
        --bedFile=article.hu_single.bed      \
        --bimFile=article.hu_single.bim      \
        --famFile=article.hu_single.fam      \
        --AlleleOrder=alt-first \
        --SAIGEOutputFile=step2 \
        --minMAF=0 \
        --minMAC=10 \
        --LOCO=FALSE \
        --GMMATmodelFile=step1.rda \
        --varianceRatioFile=step1.varianceRatio.txt \
        --is_output_moreDetails=TRUE
```

## Genetic contribution for litter size
First extract the target loci (chr6:30150129, chr6:30014787, chr6:32425024), get snp3.raw
```R
df=read.table("snp3.raw",header=T)
phe = read.table("changao.mean.txt",header = T)

colnames(df)[7:9] =c("unc5c","bmpr1b","grid2")
colnames(phe)[1]="FID"

data=merge(phe,df,by="FID",sort = F)
data=data[,-c(3,4,5,6,7)]
data$sort = ifelse(data$Mean_Value ==1,0,1)  #0 for single lamb, 1 for multiple lambs
data$sort=as.factor(data$sort)


m_0 = glm(sort ~ 1,data = data,family = "binomial")
m_logistic1=glm(sort ~ bmpr1b,data = data,family = "binomial")
m_logistic2=glm(sort ~ bmpr1b + unc5c,data = data,family = "binomial")
m_logistic3=glm(sort ~ unc5c+bmpr1b+grid2,data = data,family = "binomial")
m_logistic4=glm(sort ~ unc5c+bmpr1b+grid2+bmpr1b*grid2,data = data,family = "binomial")
m_logistic5=glm(sort ~ unc5c+bmpr1b+grid2+unc5c*grid2,data = data,family = "binomial")
m_logistic6=glm(sort ~ unc5c+bmpr1b+grid2+bmpr1b*unc5c,data = data,family = "binomial")

summary(m_logistic1)
summary(m_logistic2)
summary(m_logistic3)
summary(m_logistic4)
summary(m_logistic5)
summary(m_logistic6)

AIC(m_logistic1)
AIC(m_logistic2)
AIC(m_logistic3)
AIC(m_logistic4)
AIC(m_logistic5)
AIC(m_logistic6)

#McFadden_R2
logLik_full1 <- logLik(m_logistic1)
logLik_full2 <- logLik(m_logistic2)
logLik_full3 <- logLik(m_logistic3)
logLik_full4 <- logLik(m_logistic4)
logLik_full5 <- logLik(m_logistic5)
logLik_full6 <- logLik(m_logistic6)
logLik_null <- logLik(m_0)
1 - (logLik_full1 / logLik_null)
1 - (logLik_full2 / logLik_null)
1 - (logLik_full3 / logLik_null)
1 - (logLik_full4 / logLik_null)
1 - (logLik_full5 / logLik_null)
1 - (logLik_full6 / logLik_null)

#p value
anova(m_logistic1,m_0)
anova(m_logistic2,m_0)
anova(m_logistic3,m_0)
anova(m_logistic4,m_0)
anova(m_logistic5,m_0)
anova(m_logistic6,m_0)
```

