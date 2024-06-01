##########本程序是用来从Hista2和STAR比对后采用HTSeq计数结果一致到下游分析#####################
######运用到的文件有，
###sif文件：lcs_rnaseq20231219.sif（/public/workspace/lincs/lab7/RNAseq/downstream/中）
###hisat-count2-full文件夹下有6个计数文件
###GTF.grch38.98_geneid_genename.Rdata，该文件是通过gtf-geneid-genename.R代码运行得来，是ensemble ID和symbol对应表
###gmt文件：h.all.v7.0.symbols.gmt，h.all.v7.1.entrez.gmt
###meta文件：以后处理转录组数据必须得有meta文件
###clusterProfiler这个包有时会出bug，联网匹配数据出问题，因此有时富集不到，还有富集时不同时间结果也不同
################第一部分：读入HTSeq计数文件构建矩阵
setwd("/public/workspace/stu21230110/myNGS/bigwork/7.downstream/")#改成你的合适的路径
library(readr)
#做差异分析
library(DESeq2)
library(ggplot2)
#以下为合并HTseq生成的txt文件代码，需要在同目录下建立一个excel文件，包括列：文件名，样品名，分组信息，
#这些也是测序公司一般的操作，一般送样本时需要知道样本名称及分组信息，这样将TCGA的分析思路引入进来
#TCGA的样品信息也有Meta文件，文件名，患者barcode
#excel文件读取用readxl包
library(readxl)
meta_info = read_excel("meta3.xlsx",sheet = 1)
require(data.table)
setwd(paste0(getwd(),"/3.NMDi"))
expr_df <- fread(meta_info$file_name[1])
names(expr_df) <- c("gene_id",meta_info$sample_name[1])
#内联需要dplyr包
library(dplyr)
library(tidyr)
for (i in 2:nrow(meta_info)){  
  dfnew <- fread(meta_info$file_name[i])
  names(dfnew) <- c("gene_id",meta_info$sample_name[i])
  expr_df <- inner_join(expr_df,dfnew,by="gene_id")
}
mRNAmatrix <- expr_df[1:(length(expr_df$gene_id)-5),]
setwd("../")
if(!file.exists("3.result")){dir.create("3.result")}
save(mRNAmatrix,file = "./3.result/01-mRNAmatrix_merger.RData")
write.table(mRNAmatrix,file = "./3.result/01-mRNAmatrix_merge.txt",quote = F,row.names = F)

################第二部分：写入表达矩阵
#####那表达矩阵也需要写一个文件的,这时是counts
library(DESeq2)
#以下三行变为矩阵，但是存在同名的rownames
mRNAmatrix <- mRNAmatrix %>%
  tidyr::separate(gene_id,into = c("gene_id","drop"),sep="\\.") %>%
  dplyr::select(-drop)
###下面几行代码是将gene_id转变为gene_name
load("./GTF.grch38.98_geneid_genename.Rdata")
#采用左连接，左边的全部保留
mRNAmatrix2<-left_join(x=mRNAmatrix,y=gtf_df4,by="gene_id")
#删除没能注释的行
na_row<-which(is.na(mRNAmatrix2$gene_name)==T)
if(length(na_row)==0){
  mRNAmatrix3<-mRNAmatrix2
}else{
  mRNAmatrix3<-mRNAmatrix2[-na_row,]}

dup_row<-which(duplicated(mRNAmatrix3$gene_name)==T)
mRNAmatrix4<-mRNAmatrix3[-dup_row,]
mRNAmatrix4<-as.data.frame(mRNAmatrix4)
rownames(mRNAmatrix4)<-mRNAmatrix4$gene_name
mRNAmatrix5<-mRNAmatrix4[,-c(1,8)]
mRNAmatrix6<-as.matrix(mRNAmatrix5)

#colData
colData<-meta_info[,-1]
colData$group<-as.factor(colData$group)
condition <- factor(c(rep("C",3),rep("T",3)), levels = c("C","T")) 
dds <- DESeqDataSetFromMatrix(countData=mRNAmatrix6, 
                              colData=colData, 
                              design= ~ group)
dds_counts<-assay(dds)
#去除小数点
dds_counts<-as.data.frame(dds_counts)
colnames(dds_counts)<-colData(dds)$sample_name
dds_counts$gene_id<-rownames(dds_counts)
#写入表达矩阵，行名为基因名。
#这是将所有分析的基因都记录下来
write.csv(mRNAmatrix6,file="./3.result/02-1_all_gene_expression_matrics.csv")
save.image(file = "./3.result/02-2.RData")
################第三部分：对dds对象预测预处理
#load(file = "02-2.RData")
##数据集预过滤
nrow(dds)
dds<-dds[rowSums(counts(dds))>1,]
nrow(dds)
##rlog转换
rld<-rlog(dds,blind = F)
# 计算样本之间的距离(欧几里得距离)
sampleDists <- dist(t(assay(rld)))
# 注意要用t转换assay(rld)
# 因为dist函数要求 行是不同样本 
# 列是不同的基因是
sampleDists
# 加载下面画图所需的R包
#install.packages('pheatmap')
library(pheatmap)
#install.packages("RColorBrewer")
library(RColorBrewer)

# 把sampleDists转化为matrix格式
sampleDistMatrix <- as.matrix( sampleDists )
# 查看sampleDistMatrix行名
rownames(sampleDistMatrix) 
# 设置行名为dex的处理条件和对应的细胞系
rownames(sampleDistMatrix) <- meta_info$sample_name
colnames(sampleDistMatrix) <- meta_info$sample_name
rownames(sampleDistMatrix)

# 设置heatmap颜色
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# 相关性热图绘制
pdf(file = "./3.result/03-1Euclidean.distance.pdf",width=9,height = 9)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         show_colnames = T)
dev.off()
###计算样本间距的另一种选择是皮尔森距离
#PoissonDistance函数使用原始计算矩阵（无须标准化），要求每个样本为一行，而不是一列，需要将计数矩阵转置
#install.packages("PoiClaClu")
library(PoiClaClu)
#计算皮尔森距离
poisd<-PoissonDistance(t(counts(dds)))
samplePoisDisMatrix<-as.matrix(poisd$dd)
rownames(samplePoisDisMatrix)<-meta_info$sample_name
colnames(samplePoisDisMatrix)<-meta_info$sample_name
# 相关性热图绘制
pdf(file = "./3.result/03-2Poisson.distance.pdf",width=9,height = 9)
pheatmap(samplePoisDisMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         show_colnames = T)
dev.off()

# 另外一种样本间距离可视化方法 PCA principal components analysis
pdf(file = "./3.result/03-3PCA.pdf",width=9,height = 9)
plotPCA(rld, intgroup = c("group","sample_name"))
dev.off()
#intgroup 设置条件

#MDS图
mdsData<-data.frame(cmdscale(sampleDistMatrix))
mds<-cbind(mdsData,as.data.frame(colData(rld)))
library(ggplot2)
pdf(file = "./3.result/03-4rld.MDS.pdf",width=9,height = 9)
ggplot(mds,aes(X1,X2,color=group,shape=sample_name))+geom_point(size=3)
dev.off()
save.image(file = "./3.result/03.RData")
################第四部分：差异分析
#load(file = "./3.result/03.RData")
dds = DESeq(dds)
res = results(dds, contrast=c("group","T", "C"))
res = res[order(res$pvalue),]
head(res)
summary(res)

#去除小数点
res_df<-as.data.frame(res)
res_df2<-data.frame(gene_name=rownames(res_df))
res_df5<-cbind(res_df2,res_df)
res_df5$gene_id<-rownames(res)

#写入数据
#这是将所有分析的基因都记录下来，注意基因不为0，全为0的都删了
write.csv(res_df5,file="./3.result/04-1-all_results.csv",row.names = F)
#写入找到所有的差异基因，此处用的阈值是  0.05 和 1
table(res_df5$padj<0.05)
###注意有的padj是NA
diff_gene <-subset(res_df5, padj < 0.05 & abs(log2FoldChange) > 1)
#dim(diff_gene)
#head(diff_gene)
write.csv(diff_gene,file= "./3.result/04-2htseq_DEG_all.csv")
write.table(diff_gene,file="./3.result/04-2htseq_DEG_all.txt",sep="\t",row.names=T,quote=F)
#找到所有的差异表达的上调的基因
diff_gene3 <-subset(res_df5, padj < 0.05 & log2FoldChange > 1)
#dim(diff_gene3)
#head(diff_gene3)
write.csv(diff_gene3,file= "./3.result/04-3htseq_DEG_up.csv")
write.table(diff_gene3,file="./3.result/04-3htseq_DEG_up.txt",sep="\t",row.names=T,quote=F)
#找到所有的差异表达的下调的基因
diff_gene4 <-subset(res_df5, padj < 0.05 & log2FoldChange < -1)
#dim(diff_gene4)
#head(diff_gene4)
write.csv(diff_gene4,file= "./3.result/04-4htseq_DEG_down.csv")
write.table(diff_gene4,file="./3.result/04-4htseq_DEG_down.txt",sep="\t",row.names=T,quote=F)
save.image(file="./3.result/04.RData")
################第五部分：火山图
#load(file="./3.result/04.RData")
library(ggplot2)
library(ggrepel)
#install.packages("ggthemes")
library(ggthemes)
library(EnhancedVolcano)

res6<-diff_gene
#EnhancedVolcano(res,
#    lab = rownames(res),
#   x = 'log2FoldChange',
#    y = 'pvalue')
pdf(file = "./3.result/05-1velcano_20221208.pdf",width=9,height = 9)
EnhancedVolcano(res6,
                lab = rownames(res6),
                x = 'log2FoldChange',
                y = 'padj')
dev.off()
save.image(file = "./3.result/05.RData")
################第六部分：热图
#load(file = "./3.result/05.RData")
library(pheatmap)
#有时会有多种分组方式，那就分别告诉R，此处增加一个性别组
#按照调整过的pvalue排序，从小到大
ord<-diff_gene[order(diff_gene[,6]),]
df<-mRNAmatrix6[rownames(ord),]
#此处用前200行基因画图
#实际作图时，先筛差异基因，再用差异基因画图；或者用变化大的Top几千个基因画图
for(i in 1:dim(df)[1]){
  for(k in 1:dim(df)[2]){df[i,k]=as.numeric(df[i,k])}}
df<-as.matrix(df)
Type=c(rep("C",3),rep("T",3))
names(Type)=colnames(df)
Type=as.data.frame(Type)
pdf(file="./3.result/06-1.heatmap.pdf",height=7,width=10)
pheatmap(df[,], 
         annotation=Type, 
         scale="row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(20),
         cluster_cols =F,
         show_colnames = T,
         show_rownames = T,
         fontsize = 12,
         fontsize_row=3,
         fontsize_col=10)
dev.off()
#如果热图打不开，上面热图代码再运行一次
save.image(file = "./3.result/06.RData")
################第七部分：富集分析
#load(file = "./3.result/06.RData")
pvalueFilter=0.5           
qvalueFilter=0.5           
colorSel="qvalue"
library(clusterProfiler)
library(DOSE)
library(stringr)
library(org.Hs.eg.db)
#get the ENTREZID for the next analysis
sig.gene= diff_gene
head(sig.gene)
genes<-rownames(sig.gene)
head(genes)
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]                                 
#GO
# write.table(gene,"gene.txt",quote = F,sep = "\t")
#gene2 = read.table("gene.txt",sep = "\t",header = T,check.names = 1,row.names = 1)
#gene = as.vector(gene2[,1])
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff=pvalueFilter, 
               qvalueCutoff=qvalueFilter,
               ont="all",
               readable =T)
dim(kk)
write.table(kk,file="./3.result/07-1GO.txt",sep="\t",quote=F,row.names = F)                 
pdf(file="./3.result/07-2barplot.pdf",width = 10,height = 8)
barplot=barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(barplot)
dev.off()
pdf(file="./3.result/07-3bubble.pdf",width = 10,height = 8)
bubble=dotplot(kk,showCategory = 10,split="ONTOLOGY",orderBy = "GeneRatio", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bubble)
dev.off()
save.image(file = "./3.result/07.RData")
################第八部分：KEGG富集分析
#KEGG enrichment
load(file = "./3.result/07.RData")
library(stringr)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
pvalueFilter=0.5         
qvalueFilter=0.5       
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"}
sig.gene= diff_gene
head(sig.gene)
genes<-rownames(sig.gene)
head(genes)
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezIDs)
gene=entrezIDs[entrezIDs!="NA"] 
# write.table(gene,"gene.txt",quote = F,sep = "\t")
#gene2 = read.table("gene.txt",sep = "\t",header = T,check.names = 1,row.names = 1)
#gene = as.vector(gene2[,1])        
#kegg
kk <- enrichKEGG(gene = gene,organism = "hsa",pvalueCutoff =1, qvalueCutoff =1,use_internal_data = T)  
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file="./3.result/08-1KEGG.txt",sep="\t",quote=F,row.names = F)
showNum=30
if(nrow(KEGG)<30){showNum=nrow(KEGG)}
pdf(file="./3.result/08-2barplot.pdf",width = 10,height = 7)
barplot=barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)
print(barplot)
dev.off()
pdf(file="./3.result/08-3bubble.pdf",width = 10,height = 7)
bubble=dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)
print(bubble)
dev.off()
save.image(file="./3.result/08.RData")
######第九部分GSEA
load(file="./3.result/08.RData")
library(clusterProfiler)
library(enrichplot)
library(ggridges)
library(ggnewscale)
library(ggplot2)
allDiff<-diff_gene
#获得基因列表
gene <- rownames(allDiff)
## 转换
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(logFC=allDiff$log2FoldChange,SYMBOL = rownames(allDiff))
gene_df <- merge(gene_df,gene,by="SYMBOL")

## geneList 三部曲
## 1.获取基因logFC
geneList <- gene_df$logFC
## 2.命名
names(geneList) = gene_df$ENTREZID
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)

head(geneList)
## 读入hallmarks gene set，从哪来？
hallmarks <- read.gmt("./h.all.v7.1.entrez.gmt")
# 需要网络，大家会很拥挤，但是速度很快
y <- GSEA(geneList,TERM2GENE =hallmarks,pvalueCutoff = 0.5)

### 作图
pdf(file="./3.result/09-1.pdf",width = 20,height = 20)
cnetplot(y,foldChange = geneList)
dev.off()
y2<- setReadable(y,OrgDb=org.Hs.eg.db, keyType = "ENTREZID")
pdf(file="./3.result/09-2.pdf",width = 20,height = 20)
cnetplot(y2,showCategory = 4,
         foldChange = geneList,
         colorEdge = T)
dev.off()
### 看整体分布
pdf(file="./3.result/09-3.pdf",width = 10,height = 10)
dotplot(y,showCategory=12,split=".sign")+facet_grid(~.sign)
dev.off()
### 选择需要呈现的来作图
yd <- data.frame(y)
library(enrichplot)
pdf(file="./3.result/09-4.pdf",width = 10,height = 6)
#####注意，以后选的差异通路需要改！！！！！！！
gseaplot2(y,"HALLMARK_INTERFERON_GAMMA_RESPONSE",color = "red",pvalue_table = T)
dev.off()
pdf(file="./3.result/09-5.pdf",width = 10,height = 6)
####以后y后的1需要改，因为这个只是一个通路，只能选1
gseaplot2(y,1,color = "red",pvalue_table = T)
dev.off()
pdf(file="./3.result/09-6.pdf",width = 10,height = 6)
ridgeplot(y)
dev.off()
# pdf(file="09-7.pdf",width = 10,height = 6)
# gseaplot2(y, geneSetID = 1:3)
# dev.off()

### 自己添加文字加文字
index <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"
pdf(file="./3.result/09-8.pdf",width = 10,height = 6)
gseaplot2(y,index ,color = "green")
dev.off()

anno <- yd[index , c("enrichmentScore","NES", "pvalue", "p.adjust")]
colnames(anno)[1] <- "ES"
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")
pdf(file="./3.result/09-9.pdf",width = 10,height = 6)
gseaplot2(y,index,color = "green")+
  annotate("text",0.8, 0.8, label = lab, hjust=0, vjust=0,size = 5)
dev.off()
save.image(file="./3.result/09.RData")
#########第十部分 GSVA
load(file="./3.result/09.RData")
#options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
#install.packages("msigdbr")
library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(GSEABase)
library(limma)
library(stringr)
library(ggplot2)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
genesets <- getGmt("./h.all.v7.0.symbols.gmt")

gs<-genesets
# 预览过滤后的结果
head(gs)
# 保存到文件，方便以后重复使用
rt=mRNAmatrix6
gsym.expr <- rt
head(gsym.expr)
gsva_es <- gsva(as.matrix(gsym.expr), gs)
head(gsva_es)
group_list <- data.frame(sample = colnames(gsva_es), group = c(rep("Con", 3), rep("Treat", 3)))
head(group_list)
design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_es)
design
contrast.matrix <- makeContrasts(Treat-Con, levels = design)
# 差异分析，b vs. a
fit <- lmFit(gsva_es, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)
pathway <- str_replace(row.names(x), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = x$t)
head(df)
cutoff <- 1
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))
#按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)
pdf(file="./3.result/10-gsva-mRNA.pdf",width = 25,height = 30)
ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('palegreen3', 'snow3', 'dodgerblue4'), guide = FALSE) +   
  #画2条虚线
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, #画虚线
             size = 0.3) + #线的粗细  
  #写label
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
            size = 3, #字的大小
            hjust = "outward" ) +  #字的对齐方式
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.1, label=ID, color = group),
            size = 3, hjust = "inward") +  
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +  
  xlab("") +ylab("Treated versus Control")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除y轴
dev.off()



















