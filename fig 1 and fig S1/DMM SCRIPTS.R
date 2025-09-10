# 1. 数据分割与预处理 -----------------------------------------------------------
set.seed(12345679)  # 确保结果可重复

train_data<-BQ72_BQ40_BQ50_res5_train

# 2. 构建phyloseq对象 --------------------------------------------------------
# 准备OTU表
otu_mat <- train_data %>% 
  select(-sampleID, -Group) %>% 
  t() %>%
  as.data.frame() %>%
  mutate_all(as.numeric) %>%
  as.matrix()

# 准备分类信息
tax_mat <- data.frame(Genus = rownames(otu_mat)) %>% 
  as.matrix()

# 准备样本信息
sample_df <- data.frame(sampleID = rownames(train_data)) %>%
  column_to_rownames("sampleID")

# 创建phyloseq对象
phy_obj <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(sample_df)
)

# 3. 数据转换与核心微生物筛选 ------------------------------------------------
phy_comp <- microbiome::transform(phy_obj, "compositional")
core_taxa <- core_members(phy_comp, detection = 0.001, prevalence = 0.5)
phy_core <- prune_taxa(core_taxa, phy_obj)

# 准备DMM分析数据
dmm_input <- abundances(phy_core) %>%
  t() %>%
  as.matrix() %>%
  round(., 0)

# 4. DMM模型拟合与评估 -------------------------------------------------------
set.seed(123)
dmm_models <- lapply(1:2, dmn, count = dmm_input, verbose = TRUE)

# 模型评估指标
model_metrics <- data.frame(
  Components = 1:2,
  Laplace = sapply(dmm_models, laplace),
  AIC = sapply(dmm_models, AIC),
  BIC = sapply(dmm_models, BIC)
)

# 可视化模型选择
pdf("DMM_model_selection.pdf", width = 6, height = 4)
plot(model_metrics$Components, model_metrics$Laplace, type = "b",
     xlab = "Number of Components", ylab = "Model Fit")
lines(model_metrics$Components, model_metrics$AIC, type = "b", lty = 2)
lines(model_metrics$Components, model_metrics$BIC, type = "b", lty = 3)
legend("topright", legend = c("Laplace", "AIC", "BIC"), lty = 1:3)
dev.off()

# 5. 结果提取与保存 ----------------------------------------------------------
# 选择最佳模型
best_model <- dmm_models[[which.min(model_metrics$Laplace)]]

# 群落类型分配
community_types <- data.frame(
  sampleID = rownames(dmm_input),
  CommunityType = factor(
    apply(mixture(best_model), 1, which.max),
    levels = 1:2,
    labels = paste0("CT", 1:2)
  )
)

# 群落特征可视化
pdf("Community_type_profiles.pdf", width = 6, height = 4)
for (k in 1:ncol(fitted(best_model))) {
  fitted(best_model) %>%
    as.data.frame() %>%
    rownames_to_column("Genus") %>%
    pivot_longer(-Genus, names_to = "Cluster", values_to = "Value") %>%
    filter(Cluster == paste0("V", k)) %>%
    arrange(Value) %>%
    mutate(Genus = factor(Genus, levels = unique(Genus))) %>%
    filter(abs(Value) > quantile(abs(Value), 0.8)) %>%
    ggplot(aes(x = Genus, y = Value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Community Type", k)) %>%
    print()
}
dev.off()

# 保存关键结果
save(community_types, best_model, model_metrics, 
     file = "DMM_analysis_results.Rdata")
write.table(community_types, "community_type_assignments.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)