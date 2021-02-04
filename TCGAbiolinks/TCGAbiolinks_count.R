setwd('C:\\Users\\think\\Desktop\\biolinks')
library(TCGAbiolinks)
library(DT)
library(dplyr)
library(SummarizedExperiment)

####GDCquery参数说明：
#1、project——获取TCGA中最新的不同癌种的项目号
getGDCprojects()$project_id
#2、data.category——project中有哪些数据类型
TCGAbiolinks:::getProjectSummary("TCGA-LUAD")$data_categories
#3、data.type——这个参数受到上一个参数的影响，不同的data.category,会有不同的data.type
# 如果下载表达数据，常用的设置如下：
# #下载rna-seq转录组的表达数据
# data.type = "Gene Expresion Quantification"
# #下载miRNA表达数据数据
# data.type = "miRNA Expression Quantification"
# #下载Copy Number Variation数据
# data.type = "Copy Number Segment"
#4、workflow.type——这个参数受到上两个参数的影响，不同的data.category和不同的data.type，会有不同的workflow.type；如果下载表达数据，会有三种数据，FPKM，Counts，FPKM-UQ
#5、legacy
#主要是设置TCGA数据有两不同入口可以下载，GDC Legacy Archive 和 GDC Data Portal，以下是官方的解释两种数据Legacy or Harmonized区别：大致意思为：Legacy 数据hg19和hg18为参考基因组（老数据）而且已经不再更新了，Harmonized数据以hg38为参考基因组的数据（新数据），现在一般选择Harmonized。
# Harmonized data options (legacy = FALSE)
# Legacy archive data options (legacy = TRUE)
#6、access——筛选数据是否开放，这个一般不用设置，不开放的数据也没必要了，所以都设置成：access=“open"
#7、platform——涉及到数据来源的平台，如芯片数据，甲基化数据等等平台的筛选，一般不做设置，除非要筛选特定平台的数据
#8、file.type——如果是在GDC Legacy Archive（legacy=TRUE）下载数据的时候使用；如果在GDC Data Portal，这个参数不用设置
#9、barcode——指定要下载的样品，例如：barcode =c"TCGA-14-0736-02A-01R-2005-01""TCGA-06-0211-02A-02R-2005-01"
#10、data.format
#可以设置的选项为不同格式的文件： ("VCF", "TXT", "BAM","SVS","BCR XML","BCR SSF XML", "TSV", "BCR Auxiliary XML", "BCR OMF XML", "BCR Biotab", "MAF", "BCR PPS XML", "XLSX")，通常情况下不用设置，默认就行；
#11、experimental.strategy
# 用于过滤不同的实验方法得到的数据：
# (1) Harmonized: WXS, RNA-Seq, miRNA-Seq, Genotyping Array.
# (2) Legacy: WXS, RNA-Seq, miRNA-Seq, Genotyping Array, DNA-Seq, Methylation array, Protein expression array, WXS,CGH array, VALIDATION, Gene expression array,WGS, MSI-Mono-Dinucleotide Assay, miRNA expression array, Mixed strategies, AMPLICON, Exon array, Total RNA-Seq, Capillary sequencing, Bisulfite-Seq
#12、sample.type
#对样本的类型进行过滤，例如，原发癌组织，复发癌等等；

##########################1、检索数据
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts",
                  legacy = F,
                  access = "open",
                  # platform = ,
                  # file.type = ,
                  # barcode = ,
                  # data.format = ,
                  experimental.strategy = "RNA-Seq"
                  # sample.type = ,
)
query[1:5,1:5]
## 检束结果
results<-getResults(query)
dim(results)
results[1:5,1:5]
colnames(results)

datatable(getResults(query, cols = c("data_type","cases")),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

results<-getResults(query,cols=c("cases"))
# TP代表PRIMARY SOLID TUMOR；NT代表Solid Tissue Normal（其他组织样本可参考学习文档）
dataSmTP <- TCGAquery_SampleTypes(barcode = results,
                                  typesample = "TP")
#View(dataSmTP )# 从samplesDown中筛选出375个TP样本barcodes
dataSmNT <- TCGAquery_SampleTypes(barcode = results,
                                  typesample = "NT")

###设置barcodes参数，筛选符合要求的375个肿瘤样本数据和32正常组织数据
#barcode参数：根据传入barcode进行数据过滤
queryDown <- GDCquery(project = "TCGA-LUAD", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts", 
                      legacy = F,
                      access = "open",
                      # platform = ,
                      # file.type = ,
                      # barcode = ,
                      # data.format = ,
                      experimental.strategy = "RNA-Seq"
                      # sample.type = ,
                      barcode = c(dataSmTP, dataSmNT)
                      )
#######################2、mRNA Expression数据下载
GDCdownload(queryDown, #上面queryDown的结果
            method = "api",#两种方法，api或者gdc-client
            files.per.chunk = 6)#分为多个片段下载,可以解决下载容易中断的问题。
##下载完成，前面如果已经下载了，这里不会再次下载


# GDCdownload(query, method = "api", files.per.chunk = 100)
# # method如果设置为client 需要将gdc-client软件所在的路径添加到环境变量中，参考：gdc-client下载TCGA数据；
# # query，为GDCquery查询的结果，
# # files.per.chunk = 10,设置同时下载的数量，如果网速慢建议设置的小一些，
## directory="D:/data" 数据存储的路径；


########################3、数据整理
query = queryDown
##整理数据——GDCprepare可以自动的帮我们获得基因表达数据
data <- GDCprepare(query = query, 
                   save = TRUE, 
                   save.filename = "LUAD.RData")   #存储一下，方便下载直接读取
dim(data)
data[1:5,1:5]
data_count <- assay(data)
data_count[1:5,1:5]
dim(data_count)
write.csv(data_count,file = paste("TCGA-LUAD","Counts.csv",sep = "-"))#文件保存
#临床信息
colData(data)[1:5,1:5]
#save(data_fpkm,file="LUAD_fpkm.Rdata")

# 获取亚型信息
dataSubt <- TCGAquery_subtype(tumor = "LUAD")

# 获取临床数据
dataClin <- GDCquery_clinic(project = "TCGA-LUAD","clinical") 


###########################################################################################
dataPrep1 = data
# 去除dataPrep1中的异常值，dataPrep1数据中含有肿瘤组织和正常组织的数据
##TCGAanalyze_Preprocessing()对数据进行预处理：
#使用spearman相关系数去除数据中的异常值
#TCGAanalyze_Preprocessing(object, cor.cut = 0, filename = NULL,width = 1000,                        height = 1000, datatype = names(assays(object))[1])
# 函数功能描述：
#Array Array Intensity correlation (AAIC) and correlation boxplot to define outlier 
dataPrep2 <- TCGAanalyze_Preprocessing(object = dataPrep1,
                                       cor.cut = 0.6,
                                       datatype = "HTSeq - Counts")
#将预处理后的数据dataPrep2，写入新文件“LUAD_dataPrep.csv”
write.csv(dataPrep2,file = "LUAD_dataPrep.csv",quote = FALSE)
###########################################################################################
# TCGAtumor_purity(barcodes, estimate, absolute, lump, ihc, cpe);
#使用来自5种方法的5个估计值作为阈值对TCGA样本进行过滤，这5个值是estimate, absolute, lump, ihc, cpe;
#这里设置cpe=0.6（cpe是派生的共识度量，是将所有方法的标准含量归一化后的均值纯度水平，以使它们具有相等的均值和标准差）

#筛选肿瘤纯度大于等于60%的样本数据
purityDATA <- TCGAtumor_purity(colnames(dataPrep1), 0, 0, 0, 0, 0.6)
# filtered 为被过滤的数据， pure_barcodes是我们要的肿瘤数据
Purity.LUAD<-purityDATA$pure_barcodes
normal.LUAD<-purityDATA$filtered
#获取肿瘤纯度大于60%的375个肿瘤组织样本+32个正常组织样本,共计407个样本
puried_data <-dataPrep2[,c(Purity.LUAD,normal.LUAD)]
###########################################################################################
#基因注释,需要加载“SummarizedExperiment”包，“SummarizedExperiment container”每个由数字或其他模式的类似矩阵的对象表示。行通常表示感兴趣的基因组范围和列代表样品。
library("SummarizedExperiment")
rowData(dataPrep1)  #传入数据dataPrep1必须为 
#                SummarizedExpensembl_gene_id      external_gene_name    original_ensembl_gene_id
#                   <character>        <character>
#ENSG00000000003 ENSG00000000003             TSPAN6 ENSG00000000003.13
#ENSG00000000005 ENSG00000000005               TNMD ENSG00000000005.5
#ENSG00000000419 ENSG00000000419               DPM1 ENSG00000000419.11
#ENSG00000000457 ENSG00000000457              SCYL3 ENSG00000000457.12
#ENSG00000000460 ENSG00000000460           C1orf112 ENSG00000000460.15
#将结果写入文件“puried.STAD.cancer.csv”
rownames(puried_data)<-rowData(dataPrep1)$external_gene_name
write.csv(puried_data,file = "puried.LUAD.csv",quote = FALSE)
###########################################################################################
#进行表达矩阵标准化和过滤，得到用于差异分析的表达矩阵
#TCGAanalyze_Normalization（）# 使用EDASeq软件包标准化mRNA转录本和miRNA。
#TCGAanalyze_Normalization()执行EDASeq包中的如下功能：
#1. EDASeq::newSeqExpressionSet
#2. EDASeq::withinLaneNormalization
#3. EDASeq::betweenLaneNormalization
#4. EDASeq::counts
dataNorm <- TCGAanalyze_Normalization(tabDF = puried_data,#RNAseq表达矩阵，行代表基因，列代表样本|
                                      geneInfo = geneInfo,
                                      method = "gcContent")
#|geneInfo|关于geneLength和gcContent的20531个基因的矩阵，“geneInfoHT”和“geneInfo”可选。
#method |选择标准化的方法，基于’gcContent’ 或 ’geneLength’的标准化方法可选

#将标准化后的数据再过滤，去除掉表达量较低（count较低）的基因，得到最终的数据
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", #用于过滤较低count数的基因的方法，有’quantile’, ’varFilter’, ’filter1’, ’filter2’
                                  qnt.cut =  0.25)#过滤的阈值
str(dataFilt)
write.csv(dataFilt,file = "TCGA_LUAD_final.csv",quote = FALSE) 


########################4、差异分析
####差异分析第一种
# Which samples are Primary Tumor
dataSmTP <- TCGAquery_SampleTypes(getResults(queryDown,cols="cases"),"TP") 
# which samples are solid tissue normal
dataSmNT <- TCGAquery_SampleTypes(getResults(queryDown,cols="cases"),"NT")
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmNT],
                            mat2 = dataFilt[,dataSmTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")  

####差异分析第二种
#首先读入表达矩阵文件
dataFilt_LIHC_final <- read.csv("TCGA_LUAD_final.csv", header = T,check.names = FALSE)
# 定义行名
rownames(dataFilt_LIHC_final) <- dataFilt_LIHC_final[,1]
dataFilt_LIHC_final <- dataFilt_LIHC_final[,-1]
View(dataFilt_STAD_final)
#定义样本的分组，肿瘤组和正常组
#目前只知道根据barcode14和15位数字0—9为肿瘤，大于就为正常样本
metadata <- data.frame(colnames(dataFilt_STAD_final))#View(metadata)
for (i in 1:length(metadata[,1])) {
  num <- as.numeric(substring(metadata[i,1],14,15))
  if (num %in% seq(1,9)) {metadata[i,2] <- "Tum"}
  if (num %in% seq(10,29)) {metadata[i,2] <- "Nor"}
}
colnames(metadata)[2] <- c("fenzu")
ifelse(metadata$fenzu%in%c("Tum")==TRUE,"Tum","Nor")#发现从376行开始是正常的
# 定义肿瘤样本分组
mat1 <- dataFilt_STAD_final[,1-375]
mat1 <- log(mat1+1)
# 定义正常组织样本分组
mat2 <- dataFilt_STAD_final[,376-407]
mat2 <- log(mat2+1)
#数据准备好了，开始差异表达分析
Data_DEGs <- TCGAanalyze_DEA(mat1 = mat1,#肿瘤组织的表达矩阵
                             mat2 = mat2,#正常样本的表达矩阵
                             Cond1type = "Tumor",#mat1中的样品分组信息
                             Cond2type = "Normal",#mat2中的样品分组信息
                             pipeline="limma",#用limma包还是edgeR
                             batch.factors = c("TSS"),#批处理纠正选项
                             voom = TRUE,
                             contrast.formula = "Mycontrast=Tumor-Normal")
View(Data_DEGs)#结果展示

######################################################################################
##设置logFC，挑选表达有差异的基因进行富集分析
#Data_DEGs_high_expr <- Data_DEGs[Data_DEGs$logFC >=1,]
#Genelist <- rownames(Data_DEGs_high_expr)
#ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor", Genelist)
##富集分析结果可视化##
#TCGAvisualize_EAbarplot(tf  = rownames(ansEA$ResBP),
#                        GOBPTab = ansEA$ResBP,
#                        GOCCTab = ansEA$ResCC,
#                        GOMFTab = ansEA$ResMF,
#                        PathTab = ansEA$ResPat,
#                        nRGTab = Genelist,
#                        nBar = 10, #显示条形图的数量
#                        filename = "TCGAvisualize_EAbarplot_Output.pdf")

##富集分析                          
ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",
                                RegulonList = rownames(dataDEGs))  

TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = rownames(dataDEGs),
                        nBar = 20)

## Kaplan-Meier analysis
group1 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
group2 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))

dataSurv <- TCGAanalyze_SurvivalKM(clinical_patient = dataClin,
                                   dataGE = dataFilt,
                                   Genelist = rownames(dataDEGs),
                                   Survresult = FALSE,
                                   ThreshTop = 0.67,
                                   ThreshDown = 0.33,
                                   p.cut = 0.05, group1, group2)

##	Cox-regression analysis并画网络图							  							   
require(dnet)  # to change
org.Hs.string <- dRDataLoader(RData = "org.Hs.string")
TabCoxNet <- TCGAvisualize_SurvivalCoxNET(dataClin,
                                          dataFilt, 
                                          Genelist = rownames(dataSurv),
                                          scoreConfidence = 700,
                                          org.Hs.string = org.Hs.string,
                                          titlePlot = "Case Study n.1 dnet")

