library(lidR)
library(aRchi)
library(dplyr)
library(rgl)
library(ggplot2)
library(data.table)
input_file <- "C:/Users/xipeng/OneDrive - UGent/Desktop/paper_4/01_new_tree_model/002_guaking_aspen_001.txt"
data <- fread(input_file)
colnames(data)[1:3] <- c("X", "Y", "Z")
pc <- LAS(data)
aRchi = aRchi::build_aRchi()
aRchi = aRchi::add_pointcloud(aRchi,point_cloud = pc)   #aRCHI为QSM文件，可为空，pc为点云文件
aRchi = skeletonize_pc(aRchi, D = 0.02, cl_dist = 0.02, max_d = 0.05)
# aRchi1 = smooth_skeleton(aRchi)
node_ID <- data.frame(aRchi@QSM$node_ID)
node<- unique(node_ID)
start_result <- data.frame()  # 存储树枝角度坐标
angle_result <- data.frame()  # 存储所有数据角度
number <- length(node$aRchi.QSM.node_ID)
for (i in 2:number) {
  aRchi2 <- aRchi
  aRchi2@QSM <- subset(aRchi2@QSM, node_ID == node[i,])  # 按照 node ID 提取每对树枝
  QSM <- aRchi2@QSM
  grouped_data <- QSM %>%
    group_by(axis_ID) %>%
    group_split()
  
  if (length(grouped_data) > 1) {
    for (j in 2:length(grouped_data)) {
      group3 <- grouped_data[[1]]  # 主分支点
      group4 <- grouped_data[[j]]  # 次级分支点
      
      # 检查 group3 和 group4 的行数是否大于 0
      if (nrow(group3) > 0 && nrow(group4) > 0) {
        
        a <- group3[1,] 
        b <- group4[1,]
        group <- group4
        distance <- c(b - a)  # 次级分支需要平移的距离
        group4 <- group4[-1, ] - distance  # 平移后的次级分支
        
        # 计算角度
        angle <- circular::deg(angle3d(
          o = unlist(group3[1, 1:3]),
          a = unlist(group3[nrow(group3), 4:6]),
          b = unlist(group4[nrow(group4), 4:6])
        ))
        
        # 如果角度大于 100，提取数据并进行线性拟合
        if (!is.na(angle) && angle > 100) {
          # 提取 group3 和 group 的前三列数据
          combined_data <- rbind(group3[, 1:3], group[, 1:3])
          
          # 确保 combined_data 的列是数值型
          combined_data <- as.data.frame(lapply(combined_data, as.numeric))
          
          # 创建线性模型
          lm_model <- lm(combined_data[, 3] ~ combined_data[, 1] + combined_data[, 2])
          r_squared <- summary(lm_model)$r.squared
          
          # 如果 R^2 高于 0.9，跳过本次循环
          if (r_squared > 0.9) {
            next
          }
        }
        
        # 如果角度大于 150，打印 group3 和 group 的前 6 列
        # if (!is.na(angle) && angle > 150) {
        #   cat("Angle > 150 detected. Printing group3 and group first 6 columns:\n")
        #   cat("group3:\n")
        #   print(group3[, 1:6], digits = 7)
        #   cat("group:\n")
        #   print(group[, 1:6], digits = 7)
        # }
        
        # 保存其他结果
        start_result <- rbind(start_result, a)
        angle_result <- rbind(angle_result, angle)
        result <- cbind(start_result, angle_result)
      }
    }
  }
}

result <- na.omit(result) # 删除NA值
result <- round(result, digits = 2)
colnames(result)[17] <- "Angle"
mean(result$Angle)

Q1 <- quantile(result$Angle, 0.25)  # 第一四分位数
Q3 <- quantile(result$Angle, 0.75)  # 第三四分位数
IQR <- Q3 - Q1                    # 四分位距

# 定义异常值的上下界
lower_bound <- 0
upper_bound <- Q3 + 1 * IQR

# 识别并移除异常值
result <- result[result$Angle >= lower_bound & result$Angle <= upper_bound, ]

ggplot(result, aes(x = Angle, y = ..density..)) +
  geom_histogram(binwidth = 5, fill = "green", color = "black")
mean(result$Angle)

write.table(result, file = "C:/Users/xipeng/OneDrive - UGent/Desktop/paper_4/002_guaking_aspen_018_angle.txt", row.names = FALSE, col.names = TRUE)

