setwd('C:\\Users\\think\\Desktop\\biolinks')
library(TCGAbiolinks)
library(DT)
library(dplyr)
library(SummarizedExperiment)
#建议:下载转录组层面的数据使用hg38，下载DNA层面的数据使用hg19，因为比如做SNP分析的时候很多数据库没有hg38版本的数据，都是hg19的
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

##检索数据
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - FPKM",
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
## 检索结果
results<-getResults(query)
dim(results)  #594  29
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
                      barcode = c(dataSmTP, dataSmNT))
#######################################################################################################
GDCdownload(queryDown, #上面queryDown的结果
            method = "api",#两种方法，api或者gdc-client
            files.per.chunk = 6)#分为多个片段下载,可以解决下载容易中断的问题。
##下载完成，前面如果已经下载了，这里不会再次下载
#获取Manifest文件
getManifest(queryDown,save = FALSE) 

##mRNA Expression数据下载
GDCdownload(query, method = "api", files.per.chunk = 100)
# method如果设置为client 需要将gdc-client软件所在的路径添加到环境变量中，参考：gdc-client下载TCGA数据；
# query，为GDCquery查询的结果，
# files.per.chunk = 10,设置同时下载的数量，如果网速慢建议设置的小一些，
# directory="D:/data" 数据存储的路径；

##整理数据——GDCprepare可以自动的帮我们获得基因表达数据
data <- GDCprepare(query = query, 
                   save = TRUE, 
                   save.filename = "LUAD.RData")   #存储一下，方便下载直接读取
dim(data)
data[1:5,1:5]
#针对SummarizedExperiment对象可以获得下面3个信息
#(1)表达矩阵
data_fpkm <- assay(data)
data_fpkm[1:5,1:5]
dim(data_fpkm)
#(2)特征（一般是指基因）信息的矩阵，包括特征的元数据，例如基因所在基因组范围
rowRanges(data)
#(3)临床信息
colData(data)[1:5,1:5]
#save(data_fpkm,file="LUAD_fpkm.Rdata")


###################################临床数据下载##############################
#GDC数据库有三种形式的临床数据
# indexed clinical: a refined clinical data that is created using the XML files.
# XML files: original source of the data
# BCR Biotab: tsv files parsed from XML files
# This code will get all clinical indexed data from TCGA

##clinical indexed data下载——第一种
clinical <- GDCquery_clinic(project= "TCGA-LUAD",type = "clinical")
dim(clinical)
clinical[1:4,1:4]
clinical %>%
  head %>% 
  DT::datatable(filter = 'top', 
                options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
                rownames = FALSE)
save(clinical,file="LUAD_clinical.Rdata")
write.csv(clinical, file="TCGAbiolinks-LUAD-clinical.csv")
##clinical数据下载——第二种
clinical2<-colData(data)
#write.csv(clinical2,file="TCGAbiolinks-LUAD-clinical.csv")

##BCR Biotab临床数据下载
query <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
names(clinical.BCRtab.all)

query <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab",
                  file.type = "radiation")
GDCdownload(query)
clinical.BCRtab.radiation <- GDCprepare(query)
clinical.BCRtab.all$clinical_drug_acc  %>% 
  head  %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))


##XML clinical data下载
query <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Clinical", 
                  file.type = "xml")
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")
clinical %>% 
  datatable(filter = 'top', 
            options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
            rownames = FALSE)


##获取所有临床数据  参考：https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/clinical.html#BCR_Biotab
library(data.table)
library(dplyr)
library(regexPipes)
clinical <- TCGAbiolinks:::getGDCprojects()$project_id %>% 
  regexPipes::grep("TCGA",value=T) %>% 
  sort %>% 
  plyr::alply(1,GDCquery_clinic, .progress = "text") %>% 
  rbindlist
readr::write_csv(clinical,path = paste0("all_clin_indexed.csv"))

# This code will get all clinical XML data from TCGA
getclinical <- function(proj){
  message(proj)
  while(1){
    result = tryCatch({
      query <- GDCquery(project = proj, data.category = "Clinical",file.type = "xml")
      GDCdownload(query)
      clinical <- GDCprepare_clinic(query, clinical.info = "patient")
      for(i in c("admin","radiation","follow_up","drug","new_tumor_event")){
        message(i)
        aux <- GDCprepare_clinic(query, clinical.info = i)
        if(is.null(aux) || nrow(aux) == 0) next
        # add suffix manually if it already exists
        replicated <- which(grep("bcr_patient_barcode",colnames(aux), value = T,invert = T) %in% colnames(clinical))
        colnames(aux)[replicated] <- paste0(colnames(aux)[replicated],".",i)
        if(!is.null(aux)) clinical <- merge(clinical,aux,by = "bcr_patient_barcode", all = TRUE)
      }
      readr::write_csv(clinical,path = paste0(proj,"_clinical_from_XML.csv")) # Save the clinical data into a csv file
      return(clinical)
    }, error = function(e) {
      message(paste0("Error clinical: ", proj))
    })
  }
}
clinical <- TCGAbiolinks:::getGDCprojects()$project_id %>% 
  regexPipes::grep("TCGA",value=T) %>% sort %>% 
  plyr::alply(1,getclinical, .progress = "text") %>% 
  rbindlist(fill = TRUE) %>% setDF %>% subset(!duplicated(clinical))

readr::write_csv(clinical,path = "all_clin_XML.csv")