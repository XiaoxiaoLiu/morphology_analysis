merged<-read.csv("~/work/data/lims2/neuron_recon_2/merged.csv")

summary(merged)
names(merged)

################ remove outliers, normalize data

y=data.matrix(merged[,7:27]) #21 gl  features
m=apply(y,2,mean)
std=apply(y,2,sd)
mm=t(replicate(122,m))
stdstd=t(replicate(122,std))
zscores= (y-mm)/stdstd
zscores[zscores<-3] <- -3
zscores[zscores>3] <- 3


#############   feature based clustering



heatmap.2(t(zscores), Rowv=TRUE, Colv=TRUE, col=redgreen(85), scale="column", key=T, 
          keysize=1.5,density.info="none", trace="none",cexCol=0.9,margins=c(12,15))


##################  PVALB analysis  ################
group=merged$IsPVALB

# anova analysis to sort features
data = data.frame(y = y, group = group[])
aov.out = aov(y ~ group, data)

pvalue = vector(,21);
for( i in 1:21){
  pvalue[i] = summary(aov.out)[[i]][["Pr(>F)"]][1]
}

orderP = order(pvalue)
orderedPvalue= pvalue[orderP]
par(mar=c(15,15,15,15))
plot(1:21,-log10(orderedPvalue),"s", xaxt="n",xlab="",ylab="-log10(p value)",lwd=3)
axis(side = 1, at = 1:21,labels = names(features)[orderP], las=3,cex=3)


yy=scale(zscores, center=TRUE, scale=TRUE)
newd=yy[,orderP]

unique.cre <- as.character(unique(merged$cre_line))
idx <- match(cre,unique.cre)
cre <- as.character(merged$cre_line)
mycolor1 <- rainbow(13)[idx]

binidx<-match(idx==5,c(FALSE,TRUE)
              mycolor2 <- rainbow(2)[idx]
              
              
              heatmap.2(t(newd), Rowv=FALSE, Colv=TRUE, col=redgreen(75), scale="column", key=T, 
                        ColSideColors=mycolor1, keysize=1.5,density.info="none", trace="none",cexCol=0.9,margins=c(12,15))
              