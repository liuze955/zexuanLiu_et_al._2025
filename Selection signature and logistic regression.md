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

## lasso
```R
library(plinkFile)
library(glmnet)
library(dplyr)
library(glmnet)
library(doParallel)
options(max.print=4000000)

file_geno="ld02"

geno = readBED(file_geno)
pop = read.table("ind.group.txt",header=F)
colnames(pop) = c("ind","pop","class")
geno = geno[rownames(geno) %in% pop$ind,]

pop_order = pop[match(rownames(geno),pop$ind),]

df = cbind(pop_order$class,geno)
df = data.frame(df)
colnames(df) = paste0("a",1:ncol(df))
colnames(df)[1]  ="phe"

#sample
idx = which(df$phe == 4)

train = df[-idx,]
test = df[idx,]

train$phe = as.factor(train$phe)
test$phe = as.factor(test$phe)

#lasso model
fit = glmnet(train[,-1], train$phe, family="multinomial", nlambda=50, alpha=1)
capture.output(fit, file = "glmnet.fit.csv")

png('fit.png',width = 2000, height = 1000, units = 'px',res=200)
plot(fit, xvar="lambda", label=TRUE)
dev.off()

set.seed(123)
registerDoParallel(cores=5)
cvfit = cv.glmnet(as.matrix(train[,-1]),train$phe,family="multinomial",type.measure = "class",nfolds=9,parallel=TRUE)
stopImplicitCluster()
capture.output(cvfit, file = "glmnet.cvfit.csv")

#predict
pred_result = predict(cvfit, newx=as.matrix(train[,-1]), type="class", s="lambda.1se")
write.table(pred_result,"pred_result.txt",quote=F,row.names=F)

save.image("lasso.Rdata")
```