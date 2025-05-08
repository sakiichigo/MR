#作者来源https://zhuanlan.zhihu.com/p/377356510?utm_medium=social&utm_oi=1134785072254382080
library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
#2.2 读取差异表达基因，将基因ID从GENE_SYMBOL转换为ENTREZ_ID：


#载入差异表达数据，只需基因ID(GO,KEGG,GSEA需要)和Log2FoldChange(GSEA需要)即可
info <- read.xlsx( "/Users/ZYP/Downloads/KEGG_GO/diffexp.xlsx", rowNames = F,colNames = T)

#指定富集分析的物种库
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html

#手动加入基因
genelist_data <- read.table("E:\\sacropeniaGwas\\smrResult\\checkFdr\\all_common_genes.txt", header = FALSE, sep = "\t")
genelist=genelist_data[,1]
gene <- bitr(genelist,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)


if(FALSE){
  library(org.Hs.eg.db)
  #基因ID转换#
  keytypes(org.Hs.eg.db) #查看所有可转化类型
  darkgrey <-read.csv("magenta.csv")
  symb = mapIds(x = org.Hs.eg.db,  #id转换的比对基因组（背景基因）的物种，以人为例
                keys = genelist, #将输入的gene_name列进行数据转换
                keytype = "SYMBOL", #输入数据的类型
                column = "ENTREZID")#输出数据的类型
  symb  = na.omit(symb)  #na省略entrezid_all中不是一一对应的数据情况
  symb = data.frame(symb) #将entrezid_all变成数据框格式
  
  write.csv(symb,"darkgreygenename.csv")
  getwd()
}


#gene ID转换
gene <- bitr(info$gene_symbol,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)


#转换前后的信息变化如图:


#2.3 GO分析:
GO<-enrichGO( gene$ENTREZID,#GO富集分析
                OrgDb = GO_database,
                keyType = "ENTREZID",#设定读取的gene ID类型
                ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                pvalueCutoff = 0.8,#设定p值阈值
                qvalueCutoff = 0.8,#设定q值阈值
                readable = T)
#可执行GO查看GO分析的结果，也可以保存到本地.


#2.4 KEGG分析:

KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.8,
                 qvalueCutoff = 0.8)
#可执行KEGG查看KEGG分析的结果，也可以保存到本地.


#2.5 GSEA分析:

#由于GSEA需要差异倍数的信息即Log2FoldChange，我们先要对gene转换后的ID和读入信息进行合并。


names(info) <- c('SYMBOL','Log2FoldChange','pvalue','padj')
info_merge <- merge(info,gene,by='SYMBOL')#合并转换后的基因ID和Log2FoldChange
GSEA_input <- info_merge$Log2FoldChange
names(GSEA_input) = info_merge$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)
GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 0.05)#GSEA富集分析
#可执行GSEA查看GSEA分析的结果，也可以保存到本地.


#3. GO_KEGG_GSEA可视化:

#3.1 GO/KEGG富集柱状图+点状图:

barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
barplot(KEGG,showCategory = 10,title = 'KEGG Pathway')
dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
dotplot(KEGG)

#GO barplot (BP/CC/MF)

#KEGG barplot

#GO dotplot (BP/CC/MF)

#KEGG dotplot
#3.2 富集基因与所在功能集/通路集的关联网络图：

enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE)#基因-通路关联网络图
enrichplot::cnetplot(KEGG,circular=FALSE,colorEdge = TRUE)#circluar为指定是否环化，基因过多时建议设置为FALSE

#Connection between genes and functions from GO

#Connection between genes and pathways from KEGG
#也可以以热图形式展现关联关系:

enrichplot::heatplot(GO,showCategory = 20)#基因-通路关联热图
enrichplot::heatplot(KEGG,showCategory = 20)

#此处为了便于展示热图以小样本数据处理结果为例
#3.3 富集到的功能集/通路集之间的关联网络图：

GO2 <- pairwise_termsim(GO)
KEGG2 <- pairwise_termsim(KEGG)
enrichplot::emapplot(GO2,showCategory = 20, color = "p.adjust", layout = "kk")#通路间关联网络图
enrichplot::emapplot(KEGG2,showCategory =20, color = "p.adjust", layout = "kk")

#Connection between functions from GO

#Connection between pathways from KEGG
#3.4 保存KEGG富集到的通路至本地文件并选择通路进行展示:
  
write.table(GO$ID, file = "D:/document/bioInfo/kegg/GO_IDs.txt", #将所有KEGG富集到的通路写入本地文件查看
              append = FALSE, quote = TRUE, sep = " ",
              eol = "\n", na = "NA", dec = ".", row.names = TRUE,
              col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "")
browseKEGG(KEGG,"hsa01240")#选择其中的hsa05166通路进行展示
#还可以在展示中标记二级通路，如图中我手动高亮标记了富集通路结果中信号转导的通路。


#Pathway hsa05166
#3.5 GO富集功能网络图:
GO_BP<-enrichGO( gene$ENTREZID,#GO富集分析BP模块
                   OrgDb = GO_database,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pvalueCutoff = 0.8,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.8,
                   minGSSize = 10,
                   maxGSSize = 500,
                   readable = T)
plotGOgraph(GO_BP)#GO-BP功能网络图
GO_CC<-enrichGO( gene$ENTREZID,#GO富集分析CC模块
                 OrgDb = GO_database,
                 keyType = "ENTREZID",
                 ont = "CC",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 minGSSize = 10,
                 maxGSSize = 500,
                 readable = T)
plotGOgraph(GO_CC)#GO-CC功能网络图
GO_MF<-enrichGO( gene$ENTREZID,#GO富集分析MF模块
                 OrgDb = GO_database,
                 keyType = "ENTREZID",
                 ont = "MF",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 minGSSize = 10,
                 maxGSSize = 500,
                 readable = T)
plotGOgraph(GO_MF)#GO-MF功能网络图

# 构建数据框
df <- data.frame(ID = GO_CC@result$ID, Description = GO_CC@result$Description)
#保存
write.csv(df, "D:/document/bioInfo/kegg/GO_CC_IDs.csv", row.names = FALSE)
#GO_BP

#GO_CC

#GO_MF
#3.6 GO富集弦图:
#先要将所有GO富集到的基因集所对应的类型写入本地文件从而得到BP/CC/MF各自的起始位置如我的数据里是1，2103，2410.

genedata<-data.frame(ID=info$gene_symbol,logFC=info$log2FoldChange)
write.table(GO$ONTOLOGY, file = "/Users/ZYP/Downloads/KEGG_GO/GO_ONTOLOGYs.txt", #将所有GO富集到的基因集所对应的类型写入本地文件从而得到BP/CC/MF各自的起始位置如我的数据里是1，2103，2410
            append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
GOplotIn_BP<-GO[1:10,c(2,3,7,9)] #提取GO富集BP的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_CC<-GO[2103:2112,c(2,3,7,9)]#提取GO富集CC的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_MF<-GO[2410:2419,c(2,3,7,9)]#提取GO富集MF的前10行,提取ID,Description,p.adjust,GeneID四列
GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
GOplotIn_CC$geneID <-str_replace_all(GOplotIn_CC$geneID,'/',',')
GOplotIn_MF$geneID <-str_replace_all(GOplotIn_MF$geneID,'/',',')
names(GOplotIn_BP)<-c('ID','Term','adj_pval','Genes')#修改列名,后面弦图绘制的时候需要这样的格式
names(GOplotIn_CC)<-c('ID','Term','adj_pval','Genes')
names(GOplotIn_MF)<-c('ID','Term','adj_pval','Genes')
GOplotIn_BP$Category = "BP"#分类信息
GOplotIn_CC$Category = "CC"
GOplotIn_MF$Category = "MF"
circ_BP<-GOplot::circle_dat(GOplotIn_BP,genedata) #GOplot导入数据格式整理
circ_CC<-GOplot::circle_dat(GOplotIn_CC,genedata) 
circ_MF<-GOplot::circle_dat(GOplotIn_MF,genedata) 
chord_BP<-chord_dat(data = circ_BP,genes = genedata) #生成含有选定基因的数据框
chord_CC<-chord_dat(data = circ_CC,genes = genedata) 
chord_MF<-chord_dat(data = circ_MF,genes = genedata) 
GOChord(data = chord_BP,#弦图
        title = 'GO-Biological Process',space = 0.01,#GO Term间距
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), #上下调基因颜色
        process.label = 10) #GO Term字体大小
GOChord(data = chord_CC,title = 'GO-Cellular Component',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), 
        process.label = 10) 
GOChord(data = chord_MF,title = 'GO-Mollecular Function',space = 0.01,
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), 
        process.label = 10)

#此处为了便于展示弦图以小样本数据处理结果为例
#3.7 GO富集弦表图:
GOCircle(circ_BP) #弦表图
GOCircle(circ_CC) 
GOCircle(circ_MF) 

#GO_BP

#GO_CC

#GO_MF
#3.8 GO富集系统聚类图:
chord<-chord_dat(data = circ_BP,genes = genedata) #生成含有选定基因的数据框
GOCluster(circ_BP,GOplotIn_BP$Term) #系统聚类图
chord<-chord_dat(data = circ_CC,genes = genedata)
GOCluster(circ_CC,GOplotIn_CC$Term) 
chord<-chord_dat(data = circ_MF,genes = genedata) 
GOCluster(circ_MF,GOplotIn_MF$Term) 

#GO_BP

#GO_CC

#GO_MF
#3.9 GSEA富集图:
ridgeplot(GSEA_KEGG) 
gseaplot2(GSEA_KEGG,1)
gseaplot2(GSEA_KEGG,1:30)#30是根据ridgeplot中有30个富集通路得到的

#GSEA

#GSEA on gene set 1

#GSEA on gene sets 1-30
#4. GO/KEGG/GSEA富集分析圈:
#如下图所示为富集圈图，最外侧为富集到的功能集/通路集/基因集，其颜色表示属于不同的二级分类，向内一圈为富集到各功能集/通路集/基因集的基因数目和P值，再向内一圈为上调/下调基因的比例与数目，最内一圈为该功能集/通路集/基因集在富集分析中的富集分数，内侧为图例。


#如下图所示，将上文中得到的GO/KEGG/GSEA富集分析结果和基因ID转换结果保存到本地的result_GO.txt / result_KEGG.txt / result_GSEA.txt 和result_ID文件，结合前文提到的差异表达文件的信息，通过Python提取出富集分析圈图所需要的以下信息：

#id : 富集到的GO/KEGG/GSEA的功能集/通路集/基因集的ID
#category : 对应ID所在的二级分类，GO分析中已有BP, CC, MF三种分类，KEGG和GSEA则需要自己根据通路集/基因集的注释信息进行分类
#gene_num.min : 该功能集/通路集/基因集的最小基因数，即0
#gene_num.max : 该功能集/通路集/基因集的最大基因数，即包括的基因总数
#gene_num.rich : 富集到该功能集/通路集/基因集的基因数
#-log10Pvalue : 富集分析的P值取-log10
#up.regulated : 富集到该功能集/通路集/基因集的上调基因数
#down.regulated : 富集到该功能集/通路集/基因集的下调基因数
#rich.factor : 富集系数
#所准备的文件内容如下图所示：


result_genes.txt用于GeneSYMBOL和ENTREZID的转换

GO分析结果

#KEGG和GSEA分析结果，不同的是第一列Kind信息是自己根据注释信息进行分类添加的

#差异表达文件，在前文中已进行描述
#在Python中的参考代码如下：

#同之前的脚本一样，提供用户友好型代码，直接在前端的参数栏指定参数即可。


#在R中运行以下代码绘制富集圈:
  
#需要指定的参数如下：

#path <-'/Users/ZYP/Downloads/KEGG_GO_GSEA/' #输入文件和输出图路径
#file <- 'GSEA' #输入文件名
#第一圈的ko_color <- c(rep('#F7CC13',2), rep('#954572',1), rep('#0796E0',8)) #各二级分类的颜色和数目
#第四圈的color_assign <- c('Signaling Pathways' = '#F7CC13', 'Metabolism' = '#954572', 'Else' = '#0796E0') #各二级分类的名称和颜色
#第四圈的category_legend [labels = c('Signaling Pathways', 'Metabolism', 'Else') #各二级分类的名称 / background = c('#F7CC13', '#954572', '#0796E0') #各二级分类的颜色]
#以下为GSEA富集分析圈图
path <-'/Users/ZYP/Downloads/KEGG_GO_GSEA/' #输入文件和输出图路径
file <- 'GSEA' #输入文件名
dat <- read.delim(str_c(path, file,'.txt'), sep = '\t', row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
dat$id <- factor(rownames(dat), levels = rownames(dat))

##第一圈，绘制id
pdf(str_c(path,file,'_circlize.pdf'), width = 24, height = 12)
circle_size = unit(1, 'snpc')
circos.par(gap.degree = 0.5, start.degree = 90)
plot_data <- dat[c('id', 'gene_num.min', 'gene_num.max')]
ko_color <- c(rep('#F7CC13',2), rep('#954572',1), rep('#0796E0',8))#各二级分类的颜色和数目

circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1)
circos.track(
  ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = ko_color,
  panel.fun = function(x, y) {
    ylim = get.cell.meta.data('ycenter')#ylim、xlim
    xlim = get.cell.meta.data('xcenter')
    sector.name = get.cell.meta.data('sector.index')#sector.name
    circos.axis(h = 'top', labels.cex = 0.4, labels.niceFacing = FALSE)
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)
  } )

##第二圈，绘制富集的基因和富集p值
plot_data <- dat[c('id', 'gene_num.min', 'gene_num.rich', '-log10Pvalue')]
label_data <- dat['gene_num.rich']
p_max <- round(max(dat$'-log10Pvalue')) + 1  
colorsChoice <- colorRampPalette(c('#FF906F', '#861D30'))
color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
    ylim = get.cell.meta.data('ycenter')  
    xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),1]
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)                  
  } 
)
    
##第三圈，绘制上下调基因
dat$all.regulated <- dat$up.regulated + dat$down.regulated
dat$up.proportion <- dat$up.regulated / dat$all.regulated
dat$down.proportion <- dat$down.regulated / dat$all.regulated

dat$up <- dat$up.proportion * dat$gene_num.max
plot_data_up <- dat[c('id', 'gene_num.min', 'up')]
names(plot_data_up) <- c('id', 'start', 'end')
plot_data_up$type <- 1 

dat$down <- dat$down.proportion * dat$gene_num.max + dat$up
plot_data_down <- dat[c('id', 'up', 'down')]
names(plot_data_down) <- c('id', 'start', 'end')
plot_data_down$type <- 2 

plot_data <- rbind(plot_data_up, plot_data_down)
label_data <- dat[c('up', 'down', 'up.regulated', 'down.regulated')]
color_assign <- colorRamp2(breaks = c(1, 2), col = c('red', 'blue'))

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.08, bg.border = NA, stack = TRUE, 
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...) 
    ylim = get.cell.meta.data('cell.bottom.radius') - 0.5 
    xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),3]
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  
    xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),4]
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)
  } )

##第四圈，绘制富集因子
plot_data <- dat[c('id', 'gene_num.min', 'gene_num.max', 'rich.factor')] 
label_data <- dat['category']  
color_assign <- c('Signaling Pathways' = '#F7CC13', 'Metabolism' = '#954572', 'Else' = '#0796E0')#各二级分类的名称和颜色

circos.genomicTrack(
  plot_data, ylim = c(0, 1), track.height = 0.3, bg.col = 'gray95', bg.border = NA,
  panel.fun = function(region, value, ...) {
    sector.name = get.cell.meta.data('sector.index')  #sector.name 
    circos.genomicRect(region, value, col = color_assign[label_data[sector.name,1]], border = NA, ytop.column = 1, ybottom = 0, ...) 
    circos.lines(c(0, max(region)), c(0.5, 0.5), col = 'gray', lwd = 0.3) 
  } )

category_legend <- Legend(
  labels = c('Signaling Pathways', 'Metabolism', 'Else'),#各二级分类的名称
  type = 'points', pch = NA, background = c('#F7CC13', '#954572', '#0796E0'), #各二级分类的颜色
  labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))

updown_legend <- Legend(
  labels = c('Up-regulated', 'Down-regulated'), 
  type = 'points', pch = NA, background = c('red', 'blue'), 
  labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))

pvalue_legend <- Legend(
  col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0), 
                       colorRampPalette(c('#FF906F', '#861D30'))(6)),
  legend_height = unit(3, 'cm'), labels_gp = gpar(fontsize = 8), 
  title_gp = gpar(fontsize = 9), title_position = 'topleft', title = '-Log10(Pvalue)')

lgd_list_vertical <- packLegend(category_legend, updown_legend, pvalue_legend)
grid.draw(lgd_list_vertical)

circos.clear()
dev.off()
#运行完成在路径下得到以下GSEA富集分析圈图


GSEA富集分析圈图
#以下为KEGG富集分析圈图
#path <-'/Users/ZYP/Downloads/KEGG_GO_GSEA/'输入文件和输出图路径
#file <- 'KEGG'#输入文件名
dat <- read.delim(str_c(path, file,'.txt'), sep = '\t', row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

dat$id <- factor(rownames(dat), levels = rownames(dat))

pdf(str_c(path,file,'_circlize.pdf'), width = 24, height = 12)
circle_size = unit(1, 'snpc')
circos.par(gap.degree = 0.5, start.degree = 90)
plot_data <- dat[c('id', 'gene_num.min', 'gene_num.max')] 
ko_color <- c(rep('#F7CC13',10), rep('#954572',10), rep('#0796E0',9), rep('green', 6)) #各二级分类的颜色和数目

circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1) 
circos.track(
  ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = ko_color,  
  panel.fun = function(x, y) {
    ylim = get.cell.meta.data('ycenter')  
    xlim = get.cell.meta.data('xcenter')
    sector.name = get.cell.meta.data('sector.index')  
    circos.axis(h = 'top', labels.cex = 0.4, labels.niceFacing = FALSE) 
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  
  } )

plot_data <- dat[c('id', 'gene_num.min', 'gene_num.rich', '-log10Pvalue')]  
label_data <- dat['gene_num.rich']  
p_max <- round(max(dat$'-log10Pvalue')) + 1  
colorsChoice <- colorRampPalette(c('#FF906F', '#861D30'))  
color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  
    ylim = get.cell.meta.data('ycenter')  
    xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),1]
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  
  } )

dat$all.regulated <- dat$up.regulated + dat$down.regulated
dat$up.proportion <- dat$up.regulated / dat$all.regulated
dat$down.proportion <- dat$down.regulated / dat$all.regulated

dat$up <- dat$up.proportion * dat$gene_num.max
plot_data_up <- dat[c('id', 'gene_num.min', 'up')]
names(plot_data_up) <- c('id', 'start', 'end')
plot_data_up$type <- 1  

dat$down <- dat$down.proportion * dat$gene_num.max + dat$up
plot_data_down <- dat[c('id', 'up', 'down')]
names(plot_data_down) <- c('id', 'start', 'end')
plot_data_down$type <- 2  

plot_data <- rbind(plot_data_up, plot_data_down)
label_data <- dat[c('up', 'down', 'up.regulated', 'down.regulated')]
color_assign <- colorRamp2(breaks = c(1, 2), col = c('red', 'blue'))

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  
    ylim = get.cell.meta.data('cell.bottom.radius') - 0.5 
    xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),3]
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  
    xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),4]
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  
  } )

plot_data <- dat[c('id', 'gene_num.min', 'gene_num.max', 'rich.factor')] 
label_data <- dat['category']  
color_assign <- c('Signaling Pathways' = '#F7CC13', 'Metabolism' = '#954572', 'Cancers' = '#0796E0', 'Else' = 'green')#各二级分类的名称和颜色
circos.genomicTrack(
  plot_data, ylim = c(0, 1), track.height = 0.3, bg.col = 'gray95', bg.border = NA,  
  panel.fun = function(region, value, ...) {
    sector.name = get.cell.meta.data('sector.index')  
    circos.genomicRect(region, value, col = color_assign[label_data[sector.name,1]], border = NA, ytop.column = 1, ybottom = 0, ...)  
    circos.lines(c(0, max(region)), c(0.5, 0.5), col = 'gray', lwd = 0.3)  
  } )

category_legend <- Legend(
  labels = c('Signaling Pathways', 'Metabolism', 'Cancers', 'Else'),#各二级分类的名称
  type = 'points', pch = NA, background = c('#F7CC13', '#954572', '#0796E0'), #各二级分类的颜色
  labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
updown_legend <- Legend(
  labels = c('Up-regulated', 'Down-regulated'), 
  type = 'points', pch = NA, background = c('red', 'blue'), 
  labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
pvalue_legend <- Legend(
  col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0), 
                       colorRampPalette(c('#FF906F', '#861D30'))(6)),
  legend_height = unit(3, 'cm'), labels_gp = gpar(fontsize = 8), 
  title_gp = gpar(fontsize = 9), title_position = 'topleft', title = '-Log10(Pvalue)')
lgd_list_vertical <- packLegend(category_legend, updown_legend, pvalue_legend)
grid.draw(lgd_list_vertical)
circos.clear()
dev.off()
#运行完成在路径下得到以下KEGG富集分析圈图


#KEGG富集分析圈图
#以下为GO富集分析圈图
path <-'/Users/ZYP/Downloads/KEGG_GO_GSEA/'#输入文件和输出图路径
file <- 'GO'#输入文件名
dat <- read.delim(str_c(path, file,'.txt'), sep = '\t', row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

dat$id <- factor(rownames(dat), levels = rownames(dat))

pdf(str_c(path,file,'_circlize.pdf'), width = 24, height = 12)
circle_size = unit(1, 'snpc')
circos.par(gap.degree = 0.5, start.degree = 90)
plot_data <- dat[c('id', 'gene_num.min', 'gene_num.max')] 
ko_color <- c(rep('#F7CC13',21), rep('#954572',0), rep('#0796E0',0))#各二级分类的颜色和数目

circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1) 
circos.track(
  ylim = c(0, 1), track.height = 0.05, bg.border = NA, bg.col = ko_color, 
  panel.fun = function(x, y) {
    ylim = get.cell.meta.data('ycenter')  
    xlim = get.cell.meta.data('xcenter')
    sector.name = get.cell.meta.data('sector.index') 
    circos.axis(h = 'top', labels.cex = 0.4, labels.niceFacing = FALSE) 
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE) 
  } )

plot_data <- dat[c('id', 'gene_num.min', 'gene_num.rich', '-log10Pvalue')] 
label_data <- dat['gene_num.rich']  
p_max <- round(max(dat$'-log10Pvalue')) + 1 
colorsChoice <- colorRampPalette(c('#FF906F', '#861D30'))  
color_assign <- colorRamp2(breaks = 0:p_max, col = colorsChoice(p_max + 1))

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.08, bg.border = NA, stack = TRUE,  
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
    ylim = get.cell.meta.data('ycenter') 
    xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),1]
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  
  } )

dat$all.regulated <- dat$up.regulated + dat$down.regulated
dat$up.proportion <- dat$up.regulated / dat$all.regulated
dat$down.proportion <- dat$down.regulated / dat$all.regulated

dat$up <- dat$up.proportion * dat$gene_num.max
plot_data_up <- dat[c('id', 'gene_num.min', 'up')]
names(plot_data_up) <- c('id', 'start', 'end')
plot_data_up$type <- 1  

dat$down <- dat$down.proportion * dat$gene_num.max + dat$up
plot_data_down <- dat[c('id', 'up', 'down')]
names(plot_data_down) <- c('id', 'start', 'end')
plot_data_down$type <- 2 

plot_data <- rbind(plot_data_up, plot_data_down)
label_data <- dat[c('up', 'down', 'up.regulated', 'down.regulated')]
color_assign <- colorRamp2(breaks = c(1, 2), col = c('red', 'blue'))

circos.genomicTrackPlotRegion(
  plot_data, track.height = 0.08, bg.border = NA, stack = TRUE, 
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)  
    ylim = get.cell.meta.data('cell.bottom.radius') - 0.5 
    xlim = label_data[get.cell.meta.data('sector.index'),1] / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),3]
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)  
    xlim = (label_data[get.cell.meta.data('sector.index'),2]+label_data[get.cell.meta.data('sector.index'),1]) / 2
    sector.name = label_data[get.cell.meta.data('sector.index'),4]
    circos.text(xlim, ylim, sector.name, cex = 0.4, niceFacing = FALSE)
  } )

plot_data <- dat[c('id', 'gene_num.min', 'gene_num.max', 'rich.factor')]  
label_data <- dat['category'] 
color_assign <- c('BP' = '#F7CC13', 'CC' = '#954572', 'MF' = '#0796E0')#各二级分类的名称和颜色
circos.genomicTrack(
  plot_data, ylim = c(0, 1), track.height = 0.3, bg.col = 'gray95', bg.border = NA, 
  panel.fun = function(region, value, ...) {
    sector.name = get.cell.meta.data('sector.index')  
    circos.genomicRect(region, value, col = color_assign[label_data[sector.name,1]], border = NA, ytop.column = 1, ybottom = 0, ...) 
    circos.lines(c(0, max(region)), c(0.5, 0.5), col = 'gray', lwd = 0.3) 
    
    category_legend <- Legend(
      labels = c('BP', 'CC', 'MF'),#各二级分类的名称
      type = 'points', pch = NA, background = c('#F7CC13', '#954572', '#0796E0'),#各二级分类的颜色 
      labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
    updown_legend <- Legend(
      labels = c('Up-regulated', 'Down-regulated'), 
      type = 'points', pch = NA, background = c('red', 'blue'), 
      labels_gp = gpar(fontsize = 8), grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'))
    pvalue_legend <- Legend(
      col_fun = colorRamp2(round(seq(0, p_max, length.out = 6), 0), 
                           colorRampPalette(c('#FF906F', '#861D30'))(6)),
      legend_height = unit(3, 'cm'), labels_gp = gpar(fontsize = 8), 
      title_gp = gpar(fontsize = 9), title_position = 'topleft', title = '-Log10(Pvalue)')
    lgd_list_vertical <- packLegend(category_legend, updown_legend, pvalue_legend)
    grid.draw(lgd_list_vertical)
    circos.clear()
    dev.off()
  })
        