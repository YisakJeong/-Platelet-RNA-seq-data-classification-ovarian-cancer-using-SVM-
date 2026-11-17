# 졸업 논문 

library(tibble)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)

meta_file_path <- "/nfs-data-124/hdkjh/data/TEP/RNA/metadata/PRJNA761450/metadata_final.csv"


geneexpr_file_path <- "/home/hdjys/hdjys/gene_expression.csv"

meta_data <- read.csv(meta_file_path) # stage가 NA임
gene_expression <- read.csv(geneexpr_file_path)



table(meta_data$Group)
nrow(meta_data)

#nrow(meta_data[meta_data$Group=='AC',])

# 정상과 난소암이면서 여성만 추출
meta_data <- subset(meta_data,(Group %in% c('AC','OVCAR'))& Sex =='F')


table(meta_data$Group) # AC원래 390에서 >>217 남자 제외

table(meta_data$Stage) # stage 
table(meta_data$Age) # Age 확인 

meta_data$Stage[(meta_data$Group)=='AC'] <- 'AC' # AC의 Stage를 NA에서 AC로 채움


# set index of gene expression & meta_data
rownames(gene_expression) <- gene_expression$X  # X 에 ENSG 형태  유전자 이름을 rowname으로 설정
gene_expression$X <- NULL


rownames(meta_data) <- meta_data$Sample
# meta_data$X <- NULL


# 공통된 샘플만 추출  # 공통된 샘플 총 361개 
sample <- intersect(colnames(gene_expression),meta_data$Sample)
meta_data <- meta_data[sample,]
gene_expression <- gene_expression[,sample]


ncol(gene_expression) # 둘 다 sample 수 동일
nrow(meta_data)
# 

#install.packages('caret')
#createDataPartition(meta_data)
library(caret)


table(meta_data$Group) # AC원래 404에서 >>222 남자 제외 >> 217로 줄어듦

table(meta_data$Stage) # stage 

table(meta_data$Age) # 


all(colnames(gene_expression) == rownames(meta_data)) # row 순서 맞는지 

#### PCA plot 

# Group을 factor로 (정상을 reference level로)
meta_data$Group <- factor(meta_data$Group, levels = c("AC", "OVCAR"))

# 2. DESeqDataSet 생성
dds <- DESeqDataSetFromMatrix(
  countData = gene_expression,
  colData = meta_data,
  design = ~ Group   # 필요하면 ~ Age + Group 등으로 배치효과/공변량 보정
)

# prefiltering
#smallestGroupSize <- min(table(train_col_data[, group_label]))
keep <- rowSums(counts(dds) >= 10)>=144
dds <- dds[keep, ]

# 4. 분산안정화 변환 (샘플 수가 적으면 rlog, 많으면 vst 추천)
vsd <- vst(dds, blind = TRUE)   # 그룹 정보와 무관하게 변환 (탐색적 PCA용)

# 5. PCA plot
pcaData <- plotPCA(vsd, intgroup = "Group",returnData = TRUE)
#plotPCA(vsd, intgroup = "Group")

library(ggplot2)
##### Batch Correction >> Batch 를 지우려다가 오히려 생물학적 신호를 죽일 수도 있다.
# NKI 기관에 문제가 있음?
# 차이를 지우기 어려운 부분이>> 양쪽의 sample 
# 암정상 신호를 유지하면서 차이를 줄이는게 default인데
# 데이터 구조적으로 , 심각한 문제라면 
# VuMC, NKI  #>> 틀리는 데이터가 특정 기관에 몰려있다면 
# 정상이랑 비슷한 패턴을 가지고 있으면서,

# split 나이, stage를 우선순위로,  가능하면 NKI, VUMC도 비율을 맞추는게 좋지 않을까?
pcaData <- cbind(pcaData, meta_data[rownames(pcaData), c("Institute","Age")])

# PC1과 각 변수의 관계 시각화
ggplot(pcaData, aes(PC1, PC2, color = Institute)) + geom_point()
# ggplot(pcaData, aes(PC1, PC2, color = Age)) + geom_point()


table(meta_data$Institute, meta_data$Group)


################### batch 보정
#BiocManager::install("sva")
library(sva)





############################ 나이 층화 작업 0~45, 46~60, 61~100(young, mid, old)>>> 
meta_data$Age_group <- cut(meta_data$Age,
                           breaks = c(0, 45, 60, 100),
                           labels = c("young", "mid", "old"))

table(meta_data$Age_group) # 89, 125, 147


meta_data$strata <- case_when(
  # 정상
  meta_data$Group == "AC" ~ 
    paste("AC", meta_data$Age_group, sep = "_"),
  
  # 난소암인데 stage 모름
  meta_data$Group == "OVCAR" & is.na(meta_data$Stage) ~ 
    paste("OVCAR_unknown", meta_data$Age_group, sep = "_"),
  
  # 난소암 stage 있음
  meta_data$Group == "OVCAR" & !is.na(meta_data$Stage) ~ 
    paste(meta_data$Stage, meta_data$Age_group, sep = "_")
)

################## train test split
meta_data$strata
set.seed(42)
train_idx <- createDataPartition(meta_data$strata, p = 0.8, list = FALSE)

train_data <- meta_data[train_idx, ]
test_data  <- meta_data[-train_idx, ]

train_readcount <- gene_expression[,train_idx]
test_readcount <- gene_expression[,-train_idx]

table(train_data$Group) # train AC>>175, OVCAR 121
table(test_data$Group)  # test AC 42, OVCAR 23



meta_data$split <- "test" # 일단 전체를 다 test로 둠 
meta_data$split[train_idx] <- "train" # train만 뽑아서 train으로 수정

################# train read count로 DESeq2 돌리기



train_deseq2_data <- train_readcount
train_col_data <- train_data
rownames(train_col_data) <-train_col_data$Sample

group_label <- "Group"
train_col_data[,group_label] <- factor(ifelse(train_col_data[,group_label] %in% c("AC"),"AC","OVCAR"),ordered=F)
table(train_col_data[,group_label]) # 175, 121


# AGE 가 그룹별로 다른지 t-test
t.test(meta_data$Age[meta_data$Group=="AC"], meta_data$Age[meta_data$Group=="OVCAR"]) 

# matrix 형태로 만들어야함 
cts <- as.matrix(train_deseq2_data) #;rm(deseq2_data)
all(colnames(cts)==rownames(train_col_data))

dds1 <- DESeqDataSetFromMatrix(countData = cts,colData = train_col_data,design = formula(paste("~",group_label)))
smallestGroupSize <- min(table(train_col_data[, group_label]))
prefilter <- rowSums(counts(dds1)>=10)>=smallestGroupSize
dds1 <- dds1[prefilter,]

# train data set으로만 DESeq2 
dds1 <- DESeq(dds1)


res<- results(dds1) %>% as.data.frame
resultsNames(dds1) #"Group_OVCAR_vs_AC

#res <- res %>% filter(!is.na(res$padj))
#res$abslog2foldchange <- abs(res$log2FoldChange)
resLFC <- lfcShrink(dds1,coef=resultsNames(dds1)[2],type='apeglm') %>% as.data.frame

resLFC$stat <- res$stat
resLFC <- resLFC %>% filter(!is.na(resLFC$padj))

# volcano plot
keyvals <- ifelse(
  !is.na(resLFC$padj) & resLFC$padj <0.05 & resLFC$log2FoldChange>1.0 ,"#e41a1c",
  ifelse(!is.na(resLFC$padj) & resLFC$padj < 0.05 & resLFC$log2FoldChange< -1.0, "#377eb8",
         "#999999"))

names(keyvals)[keyvals=="#e41a1c"] <- "Up-regulated"
names(keyvals)[keyvals=="#377eb8"] <- "Down-regulated"
names(keyvals)[keyvals=="#999999"] <- "None"

resLFC$sig <- names(keyvals)
resLFC$sig <- factor(resLFC$sig, levels = c("Up-regulated", "Down-regulated", "None"))
#BiocManager::install("EnhancedVolcano")
#library(EnhancedVolcano)

ggplot(resLFC, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_bw() +
  labs(x = "log2 Fold Change", y = "-log10 adjusted p-value", 
       title = "Volcano Plot: OVCAR vs AC")





# volcano plot
library(org.Hs.eg.db)
library(AnnotationDbi)


resLFC$symbol <- mapIds(org.Hs.eg.db,
                        keys = rownames(resLFC),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")



resLFC$gene <- rownames(resLFC)
top_genes <- resLFC[resLFC$sig != "None", ]
top_genes <- top_genes[order(top_genes$padj), ][1:15, ]

ggplot(resLFC, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "None" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(
    data = top_genes,
    aes(label = symbol),
    size = 3,
    color = "black",
    max.overlaps = 20,
    box.padding = 0.4,
    segment.color = "grey50"
  ) +
  theme_bw() +
  labs(x = "log2 Fold Change", y = "-log10 adjusted p-value",
       title = "Volcano Plot: OVCAR vs AC")

# up-regulated 
up_genes <- resLFC[resLFC$log2FoldChange > 1 & resLFC$padj < 0.05&!is.na(resLFC$padj),]
nrow(up_genes)
up_gene_ranked <- up_genes[order(up_genes$stat,decreasing=TRUE),]

#down-regulated genes
down_genes <- resLFC[resLFC$log2FoldChange < -1 & resLFC$padj < 0.05&!is.na(resLFC$padj),]
nrow(down_genes)
down_genes_ranked <- down_genes[order(down_genes$stat),]

# DEG selection
sig_genes <- rownames(resLFC)[which(resLFC$padj < 0.05 & abs(resLFC$log2FoldChange) > 1& !is.na(resLFC$padj))]



# DEG 발현 패턴 heatmap으로 시각화
# vsd1 <- vst(dds1,blind=FALSE)
# mat1 <- assay(vsd1)[sig_genes,]
# 
# annot_info <- data.frame(
#   group = vsd1$Group,
#   institute = vsd1$Institute,
#   row.names = colnames(vsd1)
# )

# my_colors <- colorRampPalette(c("navy", "white", "red"))(255)

# pheatmap(mat1,
#          scale = "row",                       # 유전자별 Z-score 계산 (패턴 비교의 핵심!)
#          clustering_distance_rows = "euclidean",
#          clustering_distance_cols = "euclidean",
#          clustering_method = "complete",
#          annotation_col = annot_info,
#          show_rownames = FALSE,              # 유전자가 많으면 FALSE, 상위 30~50개면 TRUE
#          show_colnames = FALSE,
#          color = my_colors,
#          main ="Gene Expression Heatmap of DEGs (AC vs OVCAR, Train Set)")



# 히트맵 그리기, (선형 모형을 사용할거면 아직 좀 많다) >> 실제 암정상이랑 나눠지는지 검증
# SVM   DESeq2 말고도 다른 통계 기법, 
# 샘플수보다 feature수가 많다(선형모델에서는 고려하기)
# overfitting만 아니라면, 성능이 우수한가, 생물학적으로 설명이 가능한가?

# 
# 히트맵 train이 아닌 전체 샘플에서 
# 분산 안정화 


# # heatmap dendrogram 
mat <- assay(vsd)[sig_genes,]
# # annotation data frame
# 
all(rownames(colData(dds))==rownames(meta_data))

colData(dds)$split <- meta_data$split

df1 <- as.data.frame(colData(dds)[,c("Group","split","Institute")])
colnames(df1) <- c("Group","split","Institute")
rownames(df1) <- colnames(vsd)
# 
# # significant gene 뽑은 것을 바탕으로 히트맵  
heatmap <- pheatmap(mat,annotation_col=df1,clustering_method="ward.D2",scale='row',
                    show_colnames = FALSE, show_rownames = FALSE, color=colorRampPalette(c("navy", "white", "firebrick3"))(100), 
                    breaks = seq(-3, 3, length.out = 100)) # scale="row"
# all(significant_gene %in% rownames(mat1)) 
# 
# 
# ##### chi-squared test   귀무가설 히트맵 클러스터링으로 나뉜 군집과 암/정상 그룹은 서로 독립적이다
# cluster <- cutree(heatmap$tree_col, k=2)
# 
# group <- df1$Group
# names(group) <- rownames(df1)
# common <- intersect(names(cluster), names(group))
# cluster2 <- cluster[common]
# group2   <- group[common]
# 
# tab <- table(cluster2,group2)
# chisq.test(tab)    # cluster 2에서 1이 암이고 2가 정상이라고 묶임, 그리고 히트맵의 오른쪽이 암 




# train test chi-squared test 귀무가설 : train/test 분할과 암/정상 그룹은 서로 독립적이다

# tab_inst <- table(colData(dds1)$split, colData(dds1)$Group)
# chisq.test(tab_inst)   # 0.497  귀무가설기각 




###### 기관, group    # 기관(Institution)과 그룹(AC/OVCAR) 간에 연관성이 없다
# tab_inst1 <- table(colData(dds1)$Institution, colData(dds1)$Group) # table에서 숫자가 없는게 많아서 fisher's exact test
# fisher.test(tab_inst1, simulate.p.value = TRUE)   # p-value 0.0004998

# UMAP으로 그려보기 

#splir에 train test에 대해서 test에 대해서 # 학습에서 특이적인 유전자를 잡아낸걸수도 있지 않나? 
# train, test bias 에, 의한 것이 아니라 오직 암정상에 의한 신호로 구분되었다. 
# 기관에 의한 batch effect 를 넣는 것 ? 원본 노문에서 Table S1 overview가 있다. 
# 기관 1-13 하나의 프로젝트이지만 하나의 batch를 가지고 있지 않다. 
# 기관에의한 편향을 보면 좋지 않겠냐? 난소암을 제공한 곳에서 정상을 제공한 곳이 없다. 같이 제공하는 기관이 소수이다.
#뽑은 신호가 뽑은 병원에 의한 신호인지 알기가 어려울 수도 있다. 적절한 통계 검증 방법이 있는지, 2by 2가 아니기 대문에 
# batch 보정을 하든지, 하나의 기관에서 난소암 정상을 제공한다면, 거기에서 뽑아서 해서 전체를 
# 

# 2. 카이제곱 검정 (샘플 수가 적다면 fisher.test 사용)



### train data 형성 
train_set <- mat1[sig_genes,]
train_set <- train_set[,rownames(train_data)]
train_set <- t(train_set) %>% as.data.frame()


# train_set에 OVCAR인지 AC인지 라벨링 필요, 그전에 샘플 순서 일치 확인
all(rownames(train_set)==rownames(train_data))

train_set$state <- ifelse(train_data$Group=='AC',0,1)

write.csv(train_set,'train_set.csv')

table(train_set$state)


### test set분리 dds train에서만 했기 때문에 test 따로 분산 안정화 작업

#vsd_train <- vst(dds, blind = TRUE)
#mat_train <- assay(vsd_train)


### 


# test에서 분산 안정화 

dds_test <- DESeqDataSetFromMatrix(countData = test_readcount,
                                   colData   = test_data,
                                   design    = ~ 1)


vsd_test <- vst(dds_test, blind = TRUE)
mat_test <- assay(vsd_test)

test_set <- mat_test[sig_genes,]
test_set <- test_set[,rownames(test_data)]
test_set <- t(test_set) %>% as.data.frame()

test_set$state <- ifelse(test_data$Group=='AC',0,1)
table(test_set$state)
write.csv(test_set,'test_set.csv')   


# GO 분석
library(clusterProfiler)
library(org.Hs.eg.db)


# GO-upregulated genes
# up gene GO
up_gene_symbol <- rownames(up_genes)
up_entrez <- bitr(up_gene_symbol,
                  fromType="ENSEMBL",
                  toType="ENTREZID",
                  OrgDb = org.Hs.eg.db
)
ego_up <- enrichGO(gene=up_entrez$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont="BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
head(ego_up)
dotplot(ego_up,showCategory=10, title="GO BP Enrichment (Up-regulated genes)")


# down-regulated genes GO
down_gene_symbol <- rownames(down_genes)

down_entrez <- bitr(down_gene_symbol,
                    fromType = "ENSEMBL",
                    toType="ENTREZID",
                    OrgDb = org.Hs.eg.db)

ego_down <- enrichGO(gene=down_entrez$ENTREZID,
                     OrgDb=org.Hs.eg.db,
                     ont="BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

dotplot(ego_down,showCategory=10,title="GO BP Enrichment (Down-regulated genes)")



# 상위 유전자 boxplot 
down_genes_ranked
top_up_genes <- c("ENSG00000103316","ENSG00000161835","ENSG00000158578",
                  "ENSG00000184205","ENSG00000118503","ENSG00000185909","ENSG00000157654",
                  "ENSG00000151575","ENSG00000179094","ENSG00000141682")
top_down_genes <- c("ENSG00000198242","ENSG00000179144","ENSG00000171858",
                    "ENSG00000177954","ENSG00000169442","ENSG00000083845",
                    "ENSG00000112306","ENSG00000109475","ENSG00000071082",
                    "ENSG00000172809")

target_genes <- c(top_up_genes,top_down_genes)

dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)
sub_norm_counts <- norm_counts[target_genes, ]
sub_log2_norm <- log2(sub_norm_counts + 1)

sub_log2_norm_t <- as.data.frame(t(sub_log2_norm))
write.csv(sub_log2_norm_t,"log2.csv")
write.csv(meta_data,"metadata.csv")