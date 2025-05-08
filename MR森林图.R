# 加载必要的包
library(readxl)
library(dplyr)
library(forestplot)

# 读取 Excel 文件，忽略列名并手动设置列名
excel_file_path <- "D:/document/bioInfo/alps-eye/result-8/select_mr.xlsx"
data <- read_excel(excel_file_path, col_names = FALSE)
colnames(data) = c("filename", "id", "pheno", "full_name")

# OR 文件路径
or_file_path <- "D:/document/bioInfo/alps-eye/result-8/or/"

# 初始化一个空的数据框来存储所有数据
all_data <- data.frame()

# 遍历 Excel 文件的每一行
for (i in 1:nrow(data)) {
  filename <- data$filename[i]
  id <- data$id[i]
  pheno <- data$pheno[i]
  
  # 构建 OR 文件的文件名
  or_file_name <- paste0(filename, " ", id, " or.csv")
  or_file_full_path <- file.path(or_file_path, or_file_name)
  
  # 检查文件是否存在
  if (file.exists(or_file_full_path)) {
    # 读取 OR 文件
    or_data <- read.csv(or_file_full_path)
    
    # 筛选 method 列包含 "Inverse" 的记录
    filtered_or_data <- or_data %>%
      filter(grepl("Inverse", method, ignore.case = TRUE))
    
    # 处理数据，加入 or_lci95 和 or_uci95
    processed_data <- filtered_or_data %>%
      mutate(
        name = paste0("Heart failure pheno",as.numeric(gsub(".*Pheno([0-9]+).*", "\\1", filename)), " on ", pheno),
        method = method,
        snp = nsnp,
        or = or,
        pval = pval,
        or_lci95 = or_lci95,
        or_uci95 = or_uci95,
        se = se
      ) %>%
      select(name, method, snp, or, pval, or_lci95, or_uci95,se)
    
    # 将处理后的数据添加到总的数据框中
    all_data <- rbind(all_data, processed_data)
  } else {
    warning(paste("文件", or_file_full_path, "不存在。"))
  }
}

# 绘制森林图
# 提取绘图所需的数据，使用 or_lci95 和 or_uci95 作为置信区间
plot_data <- all_data %>%
  select(name, method, or, pval, or_lci95, or_uci95) %>%
  rename(lower = or_lci95, upper = or_uci95)

# 准备森林图的标签，p 值用科学计数法标识
labels <- cbind(
  plot_data$name,
  plot_data$method,
  sprintf("%.2f (%.2f - %.2f)", plot_data$or, plot_data$lower, plot_data$upper),
  sprintf("%.2e", plot_data$pval)  # 使用科学计数法格式化 p 值
)

# 添加表头
col_names <- c("Name", "Method", "OR (95% CI)", "P - value")
labels <- rbind(col_names, labels)

# 构建 hrzl_lines 列表
total_rows <- nrow(labels) + 1
hrzl_lines <- vector("list", total_rows)
hrzl_lines[[1]] <- gpar(lty = 1, lwd = 2, col = "black")  # 顶部实线加粗
hrzl_lines[[2]] <- gpar(lty = 2, col = "gray")
hrzl_lines[[total_rows]] <- gpar(lty = 1, lwd = 2, col = "black")  # 底部实线加粗

# 绘制森林图
forestplot(
  labeltext = labels,
  mean = c(NA, plot_data$or),
  lower = c(NA, plot_data$lower),
  upper = c(NA, plot_data$upper),
  title = "MR forest plot",
  xlab = "OR value",
  col = fpColors(box = "royalblue", lines = "darkblue", zero = "gray"),
  grid = structure(c(1), gp = gpar(lty = 2, col = "gray")),
  boxsize = 0.2,
  lwd.ci = 2,
  xticks = seq(0, max(plot_data$upper, na.rm = TRUE), by = 0.5),
  xlog = FALSE,
  zero = 1,
  graphwidth = unit(6, "cm"),
  graph.pos = 3,
  fn.ci_norm = fpDrawCircleCI,
  colgap = unit(0.5, "cm"),
  line.margin = 0.2,
  hrzl_lines = hrzl_lines,
  is.summary = c(TRUE, rep(FALSE, nrow(plot_data))),
  txt_gp = fpTxtGp(label = gpar(cex = 0.8), ticks = gpar(cex = 0.8), title = gpar(cex = 1)),
  lwd.xaxis = 2  # 加宽 OR 值的轴线
)
