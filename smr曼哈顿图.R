# 加载必要的包
library(qqman)

# 模拟数据
set.seed(123)
n <- 1000
chr <- sample(1:22, n, replace = TRUE)
pos <- sample(1:1e6, n, replace = TRUE)
smr <- rnorm(n, 1, 0.2)  # 模拟 SMR 值

# 创建数据框
data <- data.frame(CHR = chr, BP = pos, SMR = smr)

# 添加虚拟的 SNP 列
data$SNP <- paste0("SNP_", 1:nrow(data))

# 绘制曼哈顿图
manhattan(data,
          chr = "CHR",
          bp = "BP",
          p = "SMR",
          main = "SMR 曼哈顿图",
          ylab = "SMR 值",
          col = c("blue", "red"),
          suggestiveline = FALSE,
          genomewideline = FALSE)
