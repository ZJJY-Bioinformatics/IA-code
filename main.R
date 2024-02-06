library(tidyverse)
library(mlr3verse)
library(mlr3extralearners)
library(mlr3tuning)
library(data.table)
library(multipleROC)
library(mlr3viz)
library(future)

rm(list = ls())
# read metabolite exp and group
meta = read.table("raw_data/meta//all_sample_meta.xls",header = T,row.names = 1,stringsAsFactors = F,sep= "\t")
meta_name = read_tsv("raw_data/meta/meta_name.xls")
meta_name$name = gsub(" ","_",meta_name$name)
group_data = read.table("raw_data/meta/meta_info.xls",header = T,row.names = 1,stringsAsFactors = F,sep= "\t")
group_data_zj = group_data[group_data$center == "zhujiang",]
group_data_zj %>% select(group) -> group_data_zj
group_data_bj = group_data[group_data$center == "beijing",]
group_data_bj %>% select(group) -> group_data_bj

# normalizition
meta = log10(meta+1)
# model data tidy
meta_rf_model_all = cbind(meta[rownames(group_data_zj),],group_data_zj)
# split train and test data
set.seed(123)  # set random seed

train_indices <- sample(nrow(meta_rf_model_all), nrow(meta_rf_model_all) * 0.7)
meta_rf_model = meta_rf_model_all[train_indices,]
meta_rf_model$group = factor(meta_rf_model$group, levels = unique(meta_rf_model$group))
meta_rf_model_test = meta_rf_model_all[-train_indices,]
meta_rf_model_test$group = factor(meta_rf_model_test$group, levels = unique(meta_rf_model_test$group))

meta_rf_model %>%
  pivot_longer(-group,names_to = "meta",values_to = "exp") %>%
  group_by(group,meta) %>%
  summarise(exp_sum = mean(exp)) %>%
  group_by(meta) %>%
  arrange(exp_sum) %>%
  slice_max(order_by = exp_sum,n = 1) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(group,meta) -> meta_updown

meta_rf_test = cbind(meta[rownames(group_data_bj),],group_data_bj)
meta_rf_test$group = factor(meta_rf_test$group, levels = unique(meta_rf_test$group))

meta_rf_all = cbind(meta[rownames(group_data),],group_data %>% select(group))
meta_rf_all$group = factor(meta_rf_all$group, levels = unique(meta_rf_all$group))

set.seed(123)

# Load required libraries
library(pROC)
library(randomForest)
library(caret)

meta_rf_all <- cbind(meta[rownames(group_data),], group_data %>% select(group))
meta_rf_all$group <- factor(meta_rf_all$group, levels = unique(meta_rf_all$group))
train_control <- trainControl(method = "cv", number = 10)

model <- train(group ~ ., data = meta_rf_all, method = "rf", trControl = train_control, ntree = 1000, norm.votes = FALSE)
cv_results <- model$results
print(cv_results)
predictions <- predict(model, type = "prob")
roc_results <- roc(response = meta_rf_all$group, predictor = predictions[, "UIA"])
plot(roc_results, print.thres = "best", print.auc = TRUE)

# factor importance -------
rf_data <-  as.data.table(importance(model), keep.rownames = TRUE)
colnames(rf_data)  = c("V1","V2")
as.data.frame(rf_data) %>% 
  rename(c("id" = "V1","MeanDecreaseGini" = "V2")) %>%
  left_join(meta_name,by = "id") %>%
  left_join(meta_updown, by = c("id" = "meta")) %>% 
  arrange(desc(MeanDecreaseGini)) -> rf_data

readr::write_csv(rf_data,"new_1023//randomForest_result.csv")

colnames(rf_data)[3] = "meta"

rf_data %>%
  mutate(MeanDecreaseGini = 
           case_when(
             group == "RIA" ~ MeanDecreaseGini,
             group == "UIA" ~ -MeanDecreaseGini,
           )) %>%
  arrange(desc(MeanDecreaseGini)) -> rf_data

rf_data$meta = factor(rf_data$meta,levels = rev(rf_data$meta))


ggplot() +
  geom_segment(aes(yend = rf_data$meta, 
                   x = 0,
                   y= rf_data$meta,
                   xend = rf_data$MeanDecreaseGini,
                   color = rf_data$group
                   ), alpha = 0.5,size = 1.5)+
  geom_point(data = rf_data, 
             aes(x = MeanDecreaseGini + 0.000001, 
                 y = meta, color = group),size = 5,alpha = 0.8) +
  geom_text(data = rf_data %>% filter(group == "RIA"),
            aes(x = -1,
                y = meta,
                color = group,
                label = meta),size = 8)+
  geom_text(data = rf_data %>% filter(group == "UIA"),
            aes(x = 1,
                y = meta,
                color = group,
                label = meta),hjust = 1,size = 8)+
  scale_color_manual(values = c("#E64B35","#4DBBD5")) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())+
  labs(y = "Metabolites", x = "Mean Decrease Gini Score") -> p1


ggsave(plot = p1, "new_1023//random_forest_plot.pdf",width = 10, height = 10,dpi = 300)

# select top 5 factor
rf_data %>% 
  filter(group == "RIA") %>%
  slice_max(order_by = MeanDecreaseGini,n = 5) -> rf_data_top10

library(randomForest)
library(pROC)

selected_features_list <- rf_data_top10$id
# selected_features_list = c("12_Ketone_cornerstone_cholic_acid",
#                            "Lithocholic_acid",
#                            "Glycine_deoxycholic_acid",
#                            "Glycolithocholic_acid",
#                            "Nordeoxycholic_acid")

p_list = list()
results = data.frame()
for(i in 1:length(selected_features_list)){
  selected_features  = selected_features_list[i]
  train_data <- meta_rf_model[, c(selected_features, "group")]
  test_data <- meta_rf_test[, c(selected_features, "group")]
  model <- randomForest(group ~ ., data = train_data, ntree = 1000, norm.votes = F) 
  predictions <- as.data.frame(predict(model, test_data, type = "prob"))
  predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
  predictions$observed <- test_data$group
  predictions$feature <- test_data[[selected_features]]
  
  predictions %>% arrange(predict,feature)
  roc_results <- roc(response = predictions$observed, 
                     predictor = predictions$UIA)
  roc_ci <- ci.auc(roc_results)
  
  p_list[[i]] <- multipleROC(observed~UIA,data=predictions)
  p_list[[i]]$auc
  p_list[[i]]$cutpoint
  p_list[[i]]$sens
  
  result<- data.frame(
    Feature = rf_data_top10[rf_data_top10$id == selected_features,][["meta"]],
    AUC = p_list[[i]]$auc,
    AUC_95CI = paste(round(roc_ci[1],3),round(roc_ci[3],3),sep = "-"),
    Cutoff = p_list[[i]]$cutpoint,
    Sensitivity = p_list[[i]]$sens)
    
  results = rbind(results,result)
}

plot_ROC(p_list,
         show.points = T, 
         show.eta = F, 
         show.sens = F, 
         show.AUC = T, 
         facet = F )+
  scale_color_manual(values = xbox::need_colors(5))
ggsave("new_1023//top5_factor_roc.pdf",width = 8, height = 6)
write.csv(results,"new_1023/top5_factor_roc.csv")

library(randomForest)
library(pROC)
library(dplyr)
library(ggplot2)

# Accumulate in turn
# Select the feature of top5 in the whole queue.
rf_data %>% 
  filter(group == "RIA") %>%
  slice_max(order_by = MeanDecreaseGini,n = 5) -> rf_data_top10

results <- data.frame()
p_list = list()
for (i in 1:length(rf_data_top10$id)) {
  selected_features <- rf_data_top10$id[1:i]
  
  train_data <- meta_rf_model[, c(selected_features, "group")]
  test_data <- meta_rf_test[, c(selected_features, "group")]
  
  model <- randomForest(group ~ ., data = train_data, ntree = 1000, norm.votes = FALSE) 
  predictions <- as.data.frame(predict(model, test_data, type = "prob"))
  predictions$predict <- names(predictions)[1:2][apply(predictions[, 1:2], 1, which.max)]
  predictions$observed <- test_data$group
  
  roc_results <- roc(response = predictions$observed, predictor = predictions$UIA)
  roc_ci <- ci.auc(roc_results)
  
  p_list <- multipleROC(observed ~ UIA, data = predictions)
  
  result <- data.frame(
    NumFeatures = i,
    AUC = p_list$auc,
    AUC_95CI = paste(round(roc_ci[1], 3), round(roc_ci[3], 3), sep = "-"),
    Cutoff = p_list$cutpoint,
    Sensitivity = p_list$sens
  )
  
  results <- rbind(results, result)
}

# ROC
p <- plot_ROC(p_list,
              show.points = TRUE, 
              show.eta = FALSE, 
              show.sens = FALSE, 
              show.AUC = TRUE, 
              facet = FALSE)
p <- p + labs(title = "Cumulative ROC Curve")
p
ggsave(plot = p,"new_1023//cumulative_roc.pdf", width = 8, height = 6)
dev.off()
# -------------
write.csv(results, "new_1023//cumulative_roc_results.csv", row.names = FALSE)
ggplot(results, aes(x = NumFeatures, y = AUC))+
  geom_point()+
  ggalt::geom_xspline()+
  geom_hline(yintercept = 0.99,linetype = "dashed",color = "red")+
  annotate("text", x = 2, y = 0.992, colour = "red", size=4,label = "AUC:0.99")+
  theme_classic() -> p2
ggsave(plot = p2,"new_1023//cumulative_roc_line_plot.pdf", width = 8, height = 6)
dev.off()
dev.off()
dev.off()


# Indophenol sulfate, 2- hydroxybutyric acid and 5- hydroxytryptophan are selected.
require(randomForest)
library(multipleROC)

# selected_features = rf_data[rf_data$meta %in% c("12_Ketone_cornerstone_cholic_acid",
#                                        "Lithocholic_acid",
#                                        "Glycine_deoxycholic_acid",
#                                        "Glycolithocholic_acid",
#                                        "Nordeoxycholic_acid"),]$id

selected_features <- rf_data_top10$id[c(1,2,4)]

train_data <- meta_rf_model[, c(selected_features, "group")]
test_data <- meta_rf_test[, c(selected_features, "group")]

model <- randomForest(group ~ ., data = train_data, ntree = 1000, norm.votes = FALSE) 
predictions <- as.data.frame(predict(model, test_data, type = "prob"))
predictions$predict <- names(predictions)[1:2][apply(predictions[, 1:2], 1, which.max)]
predictions$observed <- test_data$group

roc_results <- roc(response = predictions$observed, predictor = predictions$UIA)
roc_ci <- ci.auc(roc_results)

p_list <- multipleROC(observed ~ UIA, data = predictions)
ggsave("未破裂中五个.pdf",width = 8, height = 6)

# 2023-10-20
# Supplement one combination of indophenol sulfate, canine uric acid and 5- hydroxytryptophan-----.
set.seed(123) 

sample_size <- round(0.7 * nrow(train_data))
uid = sample(1:nrow(train_data), sample_size, replace = FALSE)

library(randomForest)
library(pROC)

selected_features_list <- rf_data_top10[c(1,2,5),]$id

p_list = list()
results = data.frame()

for(i in 1:(length(selected_features_list)+1)){
  if(i != (length(selected_features_list)+1)){
    selected_features  = selected_features_list[i]
  } else {
    selected_features = selected_features_list
  }
  
  train_data <- meta_rf_model[, c(selected_features, "group")]
  test_data <- meta_rf_test[, c(selected_features, "group")] 
  
  model <- randomForest(group ~ ., data = train_data[uid,], ntree = 1000, norm.votes = F) 
  predictions <- as.data.frame(predict(model, test_data, type = "prob"))
  predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
  predictions$observed <- test_data$group
  roc_results <- roc(response = predictions$observed, 
                     predictor = predictions$UIA)
  roc_ci <- ci.auc(roc_results)
  
  p_list[[i]] <- multipleROC(observed~UIA,data=predictions)
  
  result<- data.frame(
    Feature = i,
    AUC = p_list[[i]]$auc,
    AUC_95CI = paste(round(roc_ci[1],3),round(roc_ci[3],3),sep = "-"),
    Cutoff = p_list[[i]]$cutpoint,
    Sensitivity = p_list[[i]]$sens)
  
  results = rbind(results,result)
}

plot_ROC(p_list,
         show.points = T, 
         show.eta = F, 
         show.sens = F, 
         show.AUC = T, 
         facet = F )+
  scale_color_manual(values = xbox::need_colors(5))
ggsave("s3_merge_factor_roc_test.pdf",width = 8, height = 6)

results$Feature = c(selected_features_list,"all")
write.csv(results,"s3_merge_factor_roc_test.csv")

selected_features_list <- rf_data_top10[c(1,2,5),]$id

p_list = list()
results = data.frame()

for(i in 1:(length(selected_features_list)+1)){
  if(i != (length(selected_features_list)+1)){
    selected_features  = selected_features_list[i]
  } else {
    selected_features = selected_features_list
  }
  
  train_data <- meta_rf_model[, c(selected_features, "group")]
  
  model <- randomForest(group ~ ., data = train_data[uid,], ntree = 1000, norm.votes = F) 
  predictions <- as.data.frame(predict(model, train_data[-uid,], type = "prob"))
  
  predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
  predictions$observed <- train_data[-uid,]$group
  roc_results <- roc(response = predictions$observed, 
                     predictor = predictions$UIA)
  roc_ci <- ci.auc(roc_results)
  
  p_list[[i]] <- multipleROC(observed~UIA,data=predictions)
  
  result<- data.frame(
    Feature = i,
    AUC = p_list[[i]]$auc,
    AUC_95CI = paste(round(roc_ci[1],3),round(roc_ci[3],3),sep = "-"),
    Cutoff = p_list[[i]]$cutpoint,
    Sensitivity = p_list[[i]]$sens)
  
  results = rbind(results,result)
}

plot_ROC(p_list,
         show.points = T, 
         show.eta = F, 
         show.sens = F, 
         show.AUC = T, 
         facet = F )+
  scale_color_manual(values = xbox::need_colors(5))
ggsave("s3_merge_factor_roc_train.pdf",width = 8, height = 6)

results$Feature = c(selected_features_list,"all")
write.csv(results,"s3_merge_factor_roc_train.csv")

#-------------------------------------------
require(ggpubr)
require(microbiomeutilities)
rf_data$meta = gsub("-","_",rf_data$meta)
comps <- make_pairs(rf_data$group)
xbox::chdir("meta_diff_plot")
for(i in 1:37){
  p_data = cbind(meta[rownames(group_data),],group_data)[,c("center","group",rf_data$id[i])]
  #colnames(p_data)[3]  = as.character(rf_data$meta[i])
  ggplot(p_data,aes_string(x = "group",y = rf_data$id[i]))+
    geom_violin(aes(fill = group),alpha = 0.3,size = 0.5)+
    geom_boxplot(width = 0.2)+
    geom_jitter(aes(color = group),width = 0.1)+
    scale_fill_manual(values = c("#E64B35","#4DBBD5"))+
    scale_color_manual(values = c("#E64B35","#4DBBD5"))+
    facet_wrap("center~.")+
    labs(title = as.character(rf_data$meta[i]), y = "")+
    ggprism::theme_prism()->p
  
  comps <- make_pairs(rf_data$group)
  
  p + stat_compare_means(
    comparisons = comps,
    label = "p.signif",
    tip.length = 0.05,
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
      symbols = c("xxxx", "***", "**", "*", "ns")
    ),
    method = "wilcox.test") -> p2
  
  
  ggsave(plot = p2,paste(rf_data$meta[i],"_diff_exp_in_group.pdf",sep = ""),width = 8, height = 6)
}
xbox::chdir("../")
getwd()

# ==================
library(tidyverse)
library(mlr3verse)
library(mlr3extralearners)
library(mlr3tuning)
library(data.table)
library(multipleROC)
library(mlr3viz)
library(future)

rm(list = ls())
# diff sp
sp_diff = read.csv("raw_data/sp/meta_sp_diff_1023.csv",row.names = 1, header = T)
sp_cor = sp_diff


# diff kegg
kegg_diff =  read.csv("raw_data/kegg/meta_kegg_diff_1023.csv",row.names = 1, header = T)
kegg_cor = kegg_diff[,-c(1,2)]

# diff metab
colnames(meta) = meta_name$name
meta_cor = meta[,colnames(meta) %in% c("Indole_phenol_sulfate",
                                       "5-hydroxytryptophan",
                                       "Canine_uric_acid")]

kegg_cor %>% 
  t()  %>%
  as.data.frame() -> kegg_cor 

sp_cor %>% 
  t()  %>%
  as.data.frame() -> sp_cor

group_data = read_csv("data/group_data.csv",col_names = F)
group_data =  group_data %>% drop_na()

meta_cor[group_data$X2,] -> meta_cor
rownames(meta_cor) = group_data$X1

sp_cor = sp_cor[rownames(meta_cor),]
kegg_cor = kegg_cor[rownames(meta_cor),]

write_csv(sp_cor,"data/kegg_cor_diff.csv")
write_csv(meta_cor,"data/meta_cor_all.csv")
write_csv(sp_cor,"data/sp_cor_diff.csv")

# Correlation network
source("./spearman_cor_addin.r")
# iv
iv_fct = sp_cor
# mv
mv_fct = kegg_cor
# dv
dv_fct = meta_cor

iv_fct_name <- sub("_cor","",deparse(substitute(sp_cor)))
mv_fct_name <- sub("_cor","",deparse(substitute(kegg_cor)))
dv_fct_name <- sub("_cor","",deparse(substitute(meta_cor)))

#-----main-------#
cor_spearman(iv_fct,dv_fct,c(iv_fct_name ,dv_fct_name), "spearman") -> ccor_heatdat
ccor_heatdat[[1]] %>% 
  dplyr::select(meta,sp,Cor) %>%
  pivot_wider(meta,names_from = sp,values_from = Cor) %>%
  column_to_rownames("meta") -> plot_heatmap

ccor_heatdat[[1]] %>% 
  mutate(sign = case_when(
    Pval < 0.001 ~ "***",
    Pval < 0.01 ~ "**",
    Pval <= 0.05 ~ "*",
    Pval > 0.05 ~ ""
  )) %>%
  dplyr::select(meta,sp,sign) %>%
  pivot_wider(meta,names_from = sp,values_from = sign) %>%
  column_to_rownames("meta") -> plot_heatmap_text

pheatmap::pheatmap(plot_heatmap,
                   cellwidth = 10,
                   cellheight = 10,
                   display_numbers = plot_heatmap_text,
                   angle_col = 45,
                   filename = "UIA_sp_meta_cor_heatmap.pdf")
#---
cor_spearman(iv_fct,mv_fct,c(iv_fct_name ,mv_fct_name), "pearson") -> ccor1
# RIA用 pearson
cor_spearman(mv_fct,dv_fct,c(mv_fct_name,dv_fct_name),"pearson") -> ccor2
# RIA用 pearson
cor_spearman(iv_fct,dv_fct,c(iv_fct_name,dv_fct_name),"spearman") -> ccor3

ccor1[[1]] %>% filter(Pval < 0.05) -> ccor1_test
ccor2[[1]] %>% filter(Pval < 0.05) -> ccor2_test
ccor3[[1]] %>% filter(Pval < 0.05) -> ccor3_test

m1 = ccor1_test[,1:3]
m2 = ccor2_test[,1:3]
m3 = ccor3_test[,1:3]
colnames(m1) = c("source","target","cor")
colnames(m2) = c("source","target","cor")
colnames(m3) = c("source","target","cor")

all_cor = bind_rows(m1,
                    m2,
                    m3)

m4 = all_cor %>% filter(source == "Tryptophan metabolism [PATH:ko00380]"  | 
                     target == "Tryptophan metabolism [PATH:ko00380]")


m5 = all_cor %>% filter(source == "Indole_phenol_sulfate"  | 
                     target == "Indole_phenol_sulfate") %>%
  filter(!grepl("\\[",target))

need_cor = bind_rows(m4,m5)
write.csv(need_cor,"1027_need_cor.edge.csv")

# median analysis-------
require(tidyverse)

library(mediation)
medi_result = data.frame()
test =  cor_pair_data
test %>% unique()
test  = test %>% 
  filter(kegg == "K01667") %>%
  filter(meta == "Indole_phenol_sulfate")
#test  = test %>% filter(kegg == "Tryptophan metabolism [PATH:ko00380]")

for(i in 1:nrow(test)){
  
  mData <- rbind(t(iv_fct)[test[i,][[iv_fct_name]],],
                 t(mv_fct)[test[i,][[mv_fct_name]],],
                 t(dv_fct)[test[i,][[dv_fct_name]],])
  mData = as.data.frame(t(mData))
  
  colnames(mData)[1]=test[i,][1,1]
  colnames(mData)[2]=test[i,][1,2]
  colnames(mData)[3]=test[i,][1,5]
  colnames(mData) = gsub(" ","_",colnames(mData))
  colnames(mData) = gsub("^\\d-","",colnames(mData))
  colnames(mData) = make.names(colnames(mData))
  
  id = colnames(mData)
  
  iv = id[1]
  mv = id[2]
  dv = id[3]
  
  
  fitM = lm(as.formula(paste(mv,"~",iv)), data = mData)
  
  fitY <- lm(as.formula(paste(dv,"~",iv,"+",mv)), data=mData)
  
  # fitMed <- mediate(fitM, fitY, treat=iv, mediator=mv)
  # summary(fitMedBoot)-> summ_fit
  
  fitMedBoot <- mediate(fitM, fitY, boot=TRUE, sims=999, treat=iv, mediator=mv)
  summary(fitMedBoot)-> summ_fit
  
  medi_result[i,1] = summ_fit[["d.avg.p"]]
  medi_result[i,2] = summ_fit[["d.avg"]] 
  medi_result[i,3] = summ_fit[["d.avg.ci"]][1]
  medi_result[i,4] = summ_fit[["d.avg.ci"]][2]
  medi_result[i,5] = summ_fit[["z.avg.p"]]
  medi_result[i,6] = summ_fit[["z.avg"]]
  medi_result[i,7] = summ_fit[["z.avg.ci"]][1]
  medi_result[i,8] = summ_fit[["z.avg.ci"]][2]
  medi_result[i,9] = iv
  medi_result[i,10] = mv
  medi_result[i,11] = dv
  medi_result[i,12] = test[i,][1,3]
  medi_result[i,13] = test[i,][1,4]
  medi_result[i,14] = test[i,][1,6]
  medi_result[i,15] = test[i,][1,7]
}

colnames(medi_result) = c("ACME_p",
                          "ACME_beta",
                          "ACME_beta_CI_lower",
                          "ACME_beta_CI_upper",
                          "ADE_p","ADE_beta",
                          "ADE_beta_CI_lower",
                          "ADE_beta_CI_upper",
                          "iv",
                          "mv",
                          "dv",
                          "cor_iv_mv",
                          "p_iv_mv",
                          "cor_mv_dv",
                          "p_mv_dv")

medi_result_bt = medi_result

medi_result_bt_test = medi_result_bt %>% filter(ACME_p < 0.05)

write_csv(medi_result_bt_test,"medi_result_p0.05.csv")

colnames(medi_result_bt_test)

medi_result_bt_test = medi_result_bt_test

colnames(medi_result_bt_test)

medi_result_bt_test$iv = sub("\\.t__.*","",medi_result_bt_test$iv)
medi_result_bt_test$iv = gsub(".*\\.","",medi_result_bt_test$iv)
medi_result_bt_test$mv = gsub("_\\..*","",medi_result_bt_test$mv)

write.csv(medi_result_bt_test %>% 
            # group_by(mv) %>%
            # slice_max(order_by = ACME_beta, n = 5) %>%
            ungroup(), "chonjitu2_median_test.csv")

require(ggalluvial)

ggplot(data = medi_result_bt_test %>% 
         # group_by(dv,mv) %>%
         # slice_max(order_by = ACME_beta, n = 1) %>%
         ungroup() %>%
         group_by(iv,mv,dv) %>% count(),
       aes(axis1 = iv, axis2 = mv, axis3 = dv,
           y = n)) +
  scale_x_discrete(limits = c("SP", "KEGG", "META"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = mv),show.legend = F) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  ggtitle("")

ggsave("chonjitu2_median_test.pdf", width = 10, height = 15)

write_csv(medi_result_bt_test %>% group_by(iv,mv,dv) %>%
            slice_head(n = 1),"medi_result_p0.05.csv")
