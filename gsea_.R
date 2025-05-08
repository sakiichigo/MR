
# 加载easyTCGA包
library(easyTCGA)
# 加载预处理好的炎症性肠病数据集
load(file = "G:/easyTCGA_test/gse87466.Rdata")
# 将样本分组信息转换为因子类型，并指定因子水平
group <- factor(group_list,levels = c("normal","UC"))
# 查看分组情况
table(group)
# 进行差异分析，is_count = F表示数据不是count矩阵
diff_res <- diff_analysis(exprset = exprSet, group = group, is_count = F)
# 提取limma方法的差异分析结果
diff_limma <- diff_res$deg_limma
# 删除基因名包含'/'的行（对应多个symbol的探针）
diff_limma <- diff_limma[!grepl("/",diff_limma$genesymbol),]
# 查看前6行差异分析结果
head(diff_limma)

# 加载clusterProfiler包，用于富集分析
library(clusterProfiler)
# 将基因symbol转换为ENTREZID
gene_entrezid <- bitr(geneID = diff_limma$genesymbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
# 查看前6行转换结果
head(gene_entrezid)
# 将转换后的ENTREZID与差异分析结果合并
gene_entrezid <- merge(gene_entrezid,diff_limma,by.x = "SYMBOL", by.y = "genesymbol")
# 创建基因列表，以logFC值为数值，ENTREZID为名字
genelist <- gene_entrezid$logFC
names(genelist) <- gene_entrezid$ENTREZID
# 按logFC值从大到小排序
genelist <- sort(genelist,decreasing = T)
# 查看前6个基因的logFC值
head(genelist)
#genelist
genelist_data <- read.table("E:\\sacropeniaGwas\\smrResult\\checkFdr\\all_common_genes.txt", header = FALSE, sep = "\t")
genelist=genelist_data[,1]
# 加载msigdbr包，用于从msigdb数据库获取注释集
library(msigdbr)
# 下载人类的C5注释集，选择基因集名称和ENTREZ基因ID两列
m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% dplyr::select(gs_name, entrez_gene)
# 查看前6行注释集数据
head(m_t2g)
# 进行GSEA分析
gsea_res <- GSEA(genelist, 
                 TERM2GENE = m_t2g,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 seed = 456
)
# 查看第一个富集条目中的所有基因
gsea_res[[gsea_res$ID[[1]]]]

# 加载enrichplot和ggplot2包，用于可视化
library(enrichplot)
library(ggplot2)
# 绘制峰峦图，展示核心富集基因的表达倍数变化分布
ridgeplot(gsea_res,
          showCategory = 20,
          fill = "p.adjust", 
          core_enrichment = TRUE,
          label_format = 30,
          orderBy = "NES",
          decreasing = FALSE
)+
  theme(axis.text.y = element_text(size=8))
# 选择第10到15个基因集ID
ids <- gsea_res@result$ID[10:15]
# 绘制基因集的logFC分布（密度图）
gseadist(gsea_res,
         IDs = ids,
         type="density" 
)+
  theme(legend.direction = "vertical")
# 展示基因的排序以及富集分数的变化，这里展示第1个基因集
gsearank(gsea_res,
         geneSetID = 1 
)
# 获取画图数据
aa <- gsearank(gsea_res, 1, title = gsea_res[1, "Description"],output = "table")
# 查看前6行画图数据
head(aa)
# 绘制ES图，展示第1个基因集
p <- gseaplot(gsea_res, geneSetID = 1, by = "runningScore", 
              title = gsea_res$Description[1])
# 展示图形
p
# 使用ggplot2语法修改图形标题
p+theme(plot.title = element_text(size = 8,color="red"))
# 绘制logfc标准化后排序的图形，展示第1个基因集
p <- gseaplot(gsea_res, geneSetID = 1, by = "preranked", 
              title = gsea_res$Description[1])
# 展示图形
p
# 使用ggplot2语法修改图形标题
p+theme(plot.title = element_text(size = 10,color="blue"))
# 同时绘制ES和ranked-gene-list图，展示第1个基因集
p <- gseaplot(gsea_res,geneSetID = 1,title = gsea_res$Description[1])
# 展示图形
p
# 取子集修改图形标题大小
p[[1]] <- p[[1]]+theme(plot.title = element_text(size = 6))
# 展示修改后的图形
p
# 绘制GSEA富集分析的主要可视化图形（默认3个子图），展示第1个基因集
gseaplot2(gsea_res,geneSetID = 1,title = "title",
          subplots = 1:3,
          base_size = 10)
# 绘制第1个子图，展示第1个基因集
gseaplot2(gsea_res, geneSetID = 1, subplots = 1)
# 绘制第1和第2个子图，展示第1个基因集
gseaplot2(gsea_res, geneSetID = 1, subplots = 1:2)
# 将ENTREZID转换为symbol
gsea_res_symbol <- setReadable(gsea_res,"org.Hs.eg.db","ENTREZID")
# 绘制主要可视化图形，展示第1个基因集
p <- gseaplot2(gsea_res_symbol,geneSetID = 1,
               title = gsea_res_symbol$Description[1])
# 修改第1个子图标题颜色
p[[1]] <- p[[1]]+
  theme(title = element_text(color = "red"))
# 展示修改后的图形
p
# 同时绘制第4、5、6个基因集的图形
p <- gseaplot2(gsea_res,geneSetID = 4:6)
# 展示图形
p
# 修改第1个子图的图例位置和方向
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
# 展示修改后的图形
p
# 同时绘制第4、5、6个基因集的图形，并修改颜色映射
p <- gseaplot2(gsea_res,geneSetID = 4:6,
               base_size = 10,
               color = c("#E495A5", "#86B875", "#7DB0DD")
)
# 展示图形
p
# 对第4、5、6个基因集的图形进行个性化修改
p <- gseaplot2(gsea_res,geneSetID = 4:6)
p[[1]] <- p[[1]]+scale_color_viridis_d(labels=c("lalala","heiheihei","dadada"))+
  geom_hline(yintercept = 0,color="grey75", linewidth=0.8,linetype=2)+
  theme(legend.position = "top")
p[[2]] <- p[[2]]+scale_color_viridis_d()
p[[3]] <- p[[3]]+geom_hline(yintercept = 0,color="steelblue", linewidth=0.5,linetype=2)
# 展示修改后的图形
p
# 绘制包含p值信息的图形，展示第4、5、6个基因集
p <- gseaplot2(gsea_res, geneSetID = 4:6, 
               pvalue_table = TRUE 
)
# 展示图形
p
# 提取第4、5、6个基因集的NES、P值等信息
x <- gsea_res_symbol
geneSetID <- 4:6
pd <- x[geneSetID, c( "NES","pvalue", "p.adjust")]
pd <- pd[order(rownames(pd), decreasing=FALSE),]
for (i in seq_len(ncol(pd))) {pd[, i] <- format(pd[, i], digits = 4)}
# 设置表格主题
tt <- ttheme_minimal(base_size = 10,
                     core=list(
                       fg_params=list(col=c("#F8766D","#00BA38","#619CFF"))
                     )
)
# 创建表格
tp <- tableGrob(pd,rows = NULL,theme = tt)
# 修改表格行高
tp$heights <- unit(rep(0.4,nrow(tp)),"cm") 
# 将表格添加到图形中，展示第4、5、6个基因集
p <- gseaplot2(gsea_res, geneSetID = 4:6)
p[[1]] <- p[[1]]+
  annotation_custom(tp,
                    xmin = 10000,
                    xmax = 14000,
                    ymin = 0.4,
                    ymax = 0.8
  )+
  theme(plot.title = element_text(size = 5),
        legend.position = "top",
        legend.direction = "vertical"
  )
# 展示修改后的图形
p
# 选择要展示的基因，这里从第1个基因集中随机选5个
symbol <- gsea_res_symbol[[gsea_res_symbol$ID[[1]]]]
g <- sample(symbol,5)
# 在gseaplot绘制的图形中添加基因
p <- gseaplot(gsea_res_symbol, 1, by='runningScore') 
p+geom_gsea_gene(g)
# 在gseaplot2绘制的图形中添加基因，展示第6个基因集
library(ggrepel)
p <- gseaplot2(gsea_res_symbol, geneSetID = 6)
p[[1]] <- p[[1]]+geom_gsea_gene(g, geom=geom_text_repel)+
  theme(legend.position = "top",
        legend.direction = "vertical"
  )
# 展示修改后的图形
p
# 在gseaplot2绘制的图形中添加基因，展示第4、5、6个基因集
g <- sample(symbol,3)
p <- gseaplot2(gsea_res_symbol, geneSetID = 4:6)
p[[1]] <- p[[1]]+geom_gsea_gene(mapping = aes(color= Description), 
                                g,
                                geom=geom_text_repel)+
  theme(legend.position = "top",
        legend.direction = "vertical"
  )
# 展示修改后的图形
p
# 分别从第1、2、3个基因集中随机选5个基因
g11 <- sample(gsea_res_symbol[[gsea_res_symbol$ID[1]]],5)
g22 <- sample(gsea_res_symbol[[gsea_res_symbol$ID[2]]],5)
g33 <- sample(gsea_res_symbol[[gsea_res_symbol$ID[3]]],5)
# 获取第1、2、3个基因集的描述
desc <- gsea_res_symbol$Description[1:3]
# 绘制包含不同基因集不同基因的图形，展示第1、2、3个基因集
p <- gseaplot2(gsea_res_symbol, geneSetID = 1:3)
p[[1]] <- p[[1]]  + 
  geom_gsea_gene(mapping=aes(colour = Description), g11, geom=geom_text_repel, geneSet=desc[1]) + 
  geom_gsea_gene(mapping=aes(colour = Description), g22, geom=geom_text_repel, geneSet=desc[2]) +
  geom_gsea_gene(mapping=aes(colour = Description), g33, geom=geom_text_repel, geneSet=desc[3])
# 展示修改后的图形
p