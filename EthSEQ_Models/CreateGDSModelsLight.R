library(SNPRelate)
library(data.table)
library(parallel)
lab = "UCSCCoding2"
ped = fread(paste("/scratch/sharedCO/1000Genomes/",lab,".Model.AB.ped",sep=""),sep="\t",data.table=FALSE)
map = fread(paste("/scratch/sharedCO/1000Genomes/",lab,".Model.map",sep=""),sep="\t",data.table=FALSE)
sif = fread("/elaborazioni/sharedCO/CO_Shares/Code/ethseq/Models/SIF_Sel_20170208.txt",sep="\t",data.table=FALSE)
vcf = fread(paste("/scratch/sharedCO/1000Genomes/",lab,".Model.vcf",sep=""),sep="\t",data.table=FALSE)

ped[1:10,1:10]
cat(all(vcf[,1]==map[,1]),"\n")
cat(all(vcf[,2]==map[,4]),"\n")
cat(all(vcf[,3]==map[,2]),"\n")
#res = mclapply(7:ncol(ped),function(i)
#{
#  cat(i,"\n")
#  tmp = ped[,i]
#  tmp[which(tmp=="1 1")] = "A A"
#  tmp[which(tmp=="1 2")] = "A B"
#  tmp[which(tmp=="2 2")] = "B B"
#  return(tmp)
#},mc.cores=20)

#geno = matrix(unlist(res),nrow = length(res[[1]]),byrow = F)
#ped = cbind(ped[,1:6],geno)

samples = intersect(sif$Sample,ped[,2])

ped = ped[which(ped[,2]%in%samples),]
sif = sif[which(sif$Sample%in%samples),]
all(sif$Sample==ped[,2])

samples = c(sif$Sample[sample(which(sif$Race == "EUR"),350)],
		sif$Sample[sample(which(sif$Race == "AFR"),350)],
		sif$Sample[sample(which(sif$Race == "SAS"),350)],
		sif$Sample[sample(which(sif$Race == "EAS"),350)])

#samples = readLines("/elaborazioni/sharedCO/CO_Shares/Code/ethseq/Models/LighModelSamples.txt")

cat(length(samples),"\n")

ped = ped[which(ped[,2]%in%samples),]
sif = sif[which(sif$Sample%in%samples),]
all(sif$Sample==ped[,2])

idx = order(vcf$V1,vcf$V2)
vcf = vcf[idx,]
map = map[idx,]
ped = ped[,c(1:6,idx+6)]
all(vcf$V2==map$V4)

#fwrite(ped,paste("/elaborazioni/sharedCO/CO_Shares/Code/ethseq/Models/",lab,".Model.ped",sep=""),sep="\t",quote = F,col.names = F,row.names = F,nThread=5)
#fwrite(map,paste("/elaborazioni/sharedCO/CO_Shares/Code/ethseq/Models/",lab,"l.Model.map",sep=""),sep="\t",quote = F,col.names = F,row.names = F,nThread=5)
#fwrite(vcf,paste("/elaborazioni/sharedCO/CO_Shares/Code/ethseq/VCF/",lab,".Model.vcf",sep=""),sep="\t",quote = F,col.names = F,row.names = F,nThread=5)
#fwrite(sif,paste("/elaborazioni/sharedCO/CO_Shares/Code/ethseq/Models/",lab,".Model.txt",sep=""),sep="\t",quote = F,col.names = T,row.names = F,nThread=5)

geno = ped[,7:ncol(ped)]
res = mclapply(1:ncol(geno),function(i)
{
  tmp = geno[,i]
  tmp[which(!tmp%in%c("A A","A B","B B"))] = "3"
  tmp[which(tmp=="A A")] = "0"
  tmp[which(tmp=="A B")] = "1"
  tmp[which(tmp=="B B")] = "2"
  return(as.numeric(tmp))
},mc.cores=20)

geno = matrix(as.numeric(unlist(res)),nrow = length(res[[1]]),byrow = F)

mafs = 1-apply(geno,2,function(x) (length(which(x==0))*2+length(which(x==1)))/(length(which(x!=3))*2))
idx = which(mafs>0.5)
for(i in idx)
{
  tmp = geno[,i]
  tmp[which(tmp==0)] = 4
  tmp[which(tmp==2)] = 0
  tmp[which(tmp==4)] = 2
  geno[,i] = tmp
}

snp.allele = rep("A/B",ncol(geno))
snp.allele[idx] = "B/A"

snpgdsCreateGeno(paste("/elaborazioni/sharedCO/CO_Shares/Code/ethseq/Models/",lab,".Light.Model.gds",sep=""),
                 genmat = geno,
                 sample.id = ped[,2],
                 snp.id = 1:ncol(geno),
                 snp.rs.id = map[,2],
                 snp.chromosome = map[,1],
                 snp.position = map[,4],
                 snp.allele = snp.allele,
                 snpfirstdim=FALSE)

#snpgdsPED2GDS(ped.fn = "/elaborazioni/sharedCO/CO_Shares/Code/ethseq/Models/Universal.Model.ped",
#              map.fn = "/elaborazioni/sharedCO/CO_Shares/Code/ethseq/Models/Universal.Model.map",
#              out.gdsfn = "/elaborazioni/sharedCO/CO_Shares/Code/ethseq/Models/Universal.Model.gds")

genofile <- snpgdsOpen(paste("/elaborazioni/sharedCO/CO_Shares/Code/ethseq/Models/",lab,".Light.Model.gds",sep=""),readonly = F)

# Add your sample annotation
#x = read.gdsn(index.gdsn(genofile, "sample.annot"))
sex = ped[,5]
unique(sex)
sex[which(sex==1)] = "M"
sex[which(sex==2)] = "F"
samp.annot <- data.frame(pop.group = sif$Race,sex = sex)
add.gdsn(genofile, "sample.annot", samp.annot,replace=T)

sign = paste(read.gdsn(index.gdsn(genofile, "snp.rs.id")),read.gdsn(index.gdsn(genofile, "snp.chromosome")),read.gdsn(index.gdsn(genofile, "snp.position")),sep="-")
sign.vcf = paste(vcf$V3,vcf$V1,vcf$V2,sep="-")
isort = match(sign,sign.vcf)
vcf = vcf[isort,]
sign.vcf = paste(vcf$V3,vcf$V1,vcf$V2,sep="-")
all(sign.vcf==sign)

add.gdsn(genofile, "snp.ref", vcf$V4)
add.gdsn(genofile, "snp.alt", vcf$V5)

snpgdsClose(genofile)
