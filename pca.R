#install.packages(c("ggplot2","pca3d","rgl"))

#pca analysis
setwd("D:\\R-3.6.2\\P2\\PCA")    #设置工作目录
data=read.table("783-protein.txt",header=T,sep="\t",row.names=1)   #读取表格
data=t(as.matrix(data))   #矩阵转置
data.class <- rownames(data)
data.pca <- prcomp(data, scale. = TRUE)   #PCA分析
write.table(data.pca$rotation,file="PC.xls",quote=F,sep="\t")   #输出特征向量
write.table(predict(data.pca),file="newTab.xls",quote=F,sep="\t")   #输出新表
pca.sum=summary(data.pca)
write.table(pca.sum$importance,file="importance.xls",quote=F,sep="\t")   #输出PC比重
pdf(file="pcaBarplot.pdf",width=15)   #柱状图
barplot(pca.sum$importance[2,]*100,xlab="PC",ylab="percent",col="skyblue")
dev.off()
pdf(file="pcaPlot.pdf",width=15)   #碎石图
plot(pca.sum$importance[2,]*100,type="o",col="red",xlab="PC",ylab="percent")
dev.off()

#pca 2d plot
library(ggplot2)
group=c(rep("Healthy",36),rep("Recovery",4),rep("Mild",6),rep("Severe",6))  #需要大家输入
pcaPredict=predict(data.pca)
PCA = data.frame(PCA1 = pcaPredict[,1], PCA2 = pcaPredict[,2],group=group)
PCA.mean=aggregate(PCA[,1:2],list(group=PCA$group),mean)
#定义椭圆
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }
df_ell <- data.frame()
for(g in levels(PCA$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$group==g,],
                  veganCovEllipse(cov.wt(cbind(PCA1,PCA2),
                  wt=rep(1/length(PCA1),length(PCA1)))$cov,
                  center=c(mean(PCA1),mean(PCA2))))),group=g))
}

pdf(file="PCA2d-4.pdf")
ggplot(data = PCA, aes(PCA1, PCA2)) + geom_point(size=5,aes(color = group)) +
    geom_path(data=df_ell, aes(x=PCA1, y=PCA2,colour=group), size=0.6, linetype=2)+
    #annotate("text",x=PCA.mean$PCA1,y=PCA.mean$PCA2,label=PCA.mean$group)+           #把组名标在圈上
    geom_text(label=rownames(data),color="black",size=2)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

##优化图片-样本点的标签与样本点分开
install.packages("ggrepel")
library(ggrepel)
pdf(file="PCA2d-18-无标签.pdf")
ggplot(data = PCA, aes(PCA1, PCA2)) + geom_point(size=2.6,aes(color = group)) +
  geom_path(data=df_ell, aes(x=PCA1, y=PCA2,colour=group), size=0.6, linetype=2)+
  #annotate("text",x=PCA.mean$PCA1,y=PCA.mean$PCA2,label=PCA.mean$group)+
  geom_text_repel(label=rownames(data),color="black",size=3.5,max.overlaps=50)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

#pca 3d plot
library(pca3d)
library(rgl)
pca3d(data.pca, components = 1:3, group = group, show.centroids=TRUE, show.group.labels =TRUE)  #画3d图
rgl.snapshot("pca3d.png",fmt="png")   #保存3d图形

