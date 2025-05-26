## Access to environmental data
```R
library(raster)
library(vegan)
library(factoextra)

# coordinates
site=c("XiZang","HuZhou","HeNan","JiNing",
       "QingHai","ShanDong","BinCheng","NingXia",
       "KaShi","BaYin","XinJiang","SuNiTe","TaE",
       "FuHai","HuLun")
lon=c(91.11, 120.1, 113.65, 116.59, 
      101.74, 117, 118.02, 106.27, 77.65, 
      83.82, 87.68, 113.66, 85.68, 87.49, 119.76)
lat=c(29.97, 30.86, 34.76, 35.41, 36.56, 
      36.65, 37.43, 38.47,38.9, 42.89,
      43.77, 43.85, 46.5, 47.11, 49.21)
samples = data.frame(site, lon, lat)
data = samples

# extract bio1-bio19
for (i in 1:19) {
  file_name=paste0("wc2.1_5m_bio/wc2.1_5m_bio_",i,".tif")
  df = raster(file_name)
  newcol=paste0("bio",i)
  data[[newcol]]=extract(df, samples[,2:3])
}

row.names(data)=data$site
# extract elevation
elevation = raster("wc2.1_5m_elev/wc2.1_5m_elev.tif")
data$elevation= extract(elevation,samples[,2:3])

data_stander=cbind(data[,1:3],scale(data[,c(-1,-2,-3)]))

info = read.table("sample.pop.site")
data_pop_env=merge(info,data_stander,by=c("site"),all.x=T,sort=F)

# Covariance screening
library(usdm)
env=data_pop_env[,c(-1,-2,-3)]
vifstep(env,th=10)

# balance sample
sample_equal = read.table("equal.sample.txt")
env_equal = data_pop_env[data_pop_env$ind %in% sample_equal$ind,]

write.table(data_pop_env[,c("ind","pop","site","bio7","bio10","bio13","bio14","bio15")],"ind.pop.env.txt",quote = F,row.names = F)
write.table(env_equal[,c("ind","pop","site","bio7","bio10","bio13","bio14","bio15")],"ind.popEqual.env.txt",quote = F,row.names = F)
```

## RDA

Converting genotypes to matrix files
```shell
awk '{print $1,$1}' ind.popEqual.env.txt > env.keep
plink --bfile equal --keep env.keep--recode A --out equal --sheep --keep-allele-order
```
Data is large, pre-read in
```R
library(data.table)

env_path="ind.popEqual.env.txt"
env=read.table(env_path,header=T) 
colnames(env)[1] = "FID"

geno_path="equal.raw"
geno_data=fread(geno_path)

ind_sort_order = match(geno_data$FID,env$FID)
env_order=env[ind_sort_order,]

rname_ind=geno_data$FID
geno_data=geno_data[,c(-1,-2,-3,-4,-5,-6)]
env_data=env_order[,c(-1,-2,-3)]
row.names(geno_data)=rname_ind
row.names(env_data)=rname_ind
save.image("rda.data.Rdata")
```
calculate rda
```R
library(vegan)
library(data.table)

load(file="rda.data.Rdata")

# DCA
dca=decorana(geno_data)
dca1 <- max(dca$rproj[,1])
dca2 <- max(dca$rproj[,2])
dca3 <- max(dca$rproj[,3])
dca4 <- max(dca$rproj[,4])
GL <- data.frame(DCA1 = c(dca1), DCA2 = c(dca2), DCA3 = c(dca3), DCA4 = c(dca4))
write.table(GL,"out.dca",quote=F,row.names=F)

# RDA
rda <- rda(geno_data,env_data,scale = F)
RsquareAdj(rda)
```

## LFMM
```R
library(lfmm)
library(data.table)
library(qqman)

load(file="rda.data.Rdata")

env_data = env_data[,c(1,4)]
mod.lfmm = lfmm_ridge(Y = geno_data,X = env_data,K = 3)  
test.lfmm = lfmm_test(Y = geno_data,X = env_data,lfmm = mod.lfmm,calibrate = "gif")
p=test.lfmm$calibrated.pvalue
p=as.data.frame(p)
p$snp=row.names(p)

p$chr=sapply(strsplit(p$snp, "_"), function(x) x[1])
p$chr[p$chr == "X"] = 27
p$pos=sapply(strsplit(p$snp, "_"), function(x) x[2])
p$pos=as.numeric(p$pos)
p$chr=as.numeric(p$chr)
p$id=1:nrow(p)

p2=p
p2[is.na(p2)]=1
save.image("lfmm.Rdata")
fwrite(p2,"p.csv",quote=F,row.names=F)

bio=c("bio7","bio14")
#输出top位点
for (i in bio) {
name=paste0("top_",i)
a=p2[-log10(p2[,i]) > -log10(0.05/nrow(p2)) ,]
write.table(a[,c("chr","pos")],name,quote=F,row.names=F,col.names=F)
}
```

## LFMM of 90% random samples
```R
library(lfmm)
library(data.table)
library(qqman)

load("lfmm.Rdata")

#cv 10x
n = nrow(geno_data)
set.seed(123)
shuffled_indices = sample(1:n) 
folds <- cut(shuffled_indices, breaks = 10, labels = FALSE)

for (i in 1:10) {
	test_index = which(folds == i)
	train_index = which(folds != i)

	env_data_test = env_data[test_index,c(1,4)]
	geno_data_test = geno_data[test_index,]

	mod.lfmm = lfmm_ridge(Y = geno_data_test,X = env_data_test,K = 3)
	test.lfmm = lfmm_test(Y = geno_data_test,X = env_data_test,lfmm = mod.lfmm,calibrate = "gif")

	p=test.lfmm$calibrated.pvalue
	p=as.data.frame(p)
	p$snp=row.names(p)
	
	p$chr=sapply(strsplit(p$snp, "_"), function(x) x[1])
	p$chr[p$chr == "X"] = 27
	p$pos=sapply(strsplit(p$snp, "_"), function(x) x[2])
	p$pos=as.numeric(p$pos)
	p$chr=as.numeric(p$chr)
	p$id=1:nrow(p)
	

	p2=p
	cat("Number of NA values before replacement:", sum(is.na(p2)), "\n")
	p2[is.na(p2)]=1
	p2[p2==0]=1

	name_outP = paste0("lfmm.p.folds",i)
        fwrite(p2,name_outP,quote=F,row.names=F)

	#lambda
	chi_squared1 <- qchisq(1 - p2$bio7, 1)
	chi_squared2 <- qchisq(1 - p2$bio14, 1)
	lambda1 <- median(chi_squared1) / qchisq(0.5, 1)
	lambda2 <- median(chi_squared2) / qchisq(0.5, 1)
	cat("lambda_bio7",lambda1)
	cat("lambda_bio14",lambda2)

	for (bio in c("bio7","bio14")) {
		#top sites
		name_outTop = paste0("top_",bio,"_folds",i)
		a=p2[-log10(p2[,bio]) > -log10(0.05/nrow(p2)) ,]
		write.table(a[,c("chr","pos")],name_outTop,quote=F,row.names=F,col.names=F)
		#draw
		png_name1 = paste0("manhan.",bio,"folds",i,".png")
		png_name2 = paste0("qq.",bio,"folds",i,".png")
		png(png_name1,width = 2000, height = 1000, units = 'px',res=200)
		manhattan(p2 ,chr="chr", logp = T,
          			ylab="-logp",cex = 0.8, cex.axis = 0.8, suggestiveline=-log10(0.05/nrow(p2)),genomewideline=F,
          			bp="pos", p=bio ,snp = "snp" )
		dev.off()

		png(png_name2,width = 1000, height = 1000, units = 'px',res=200)
		qq(p2[,bio])
		dev.off()
	}

}
```