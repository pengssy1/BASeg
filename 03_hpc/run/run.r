
input_file ="/data/gent/vo/000/gvo00074/Xipeng/creat/Inputs/guaking_aspen_003_skeleton.txt"
output_file="/data/gent/vo/000/gvo00074/Xipeng/creat/Outputs/guaking_aspen_003_skeleton.txt.txt"

library(data.table)
library(Rfast)
library(FNN)
library(fastcluster)
library(dplyr)


# input_file = "/user/gent/447/vsc44746/LIDAR/tree1-10m-002-4scans_point.las"
# output_file = "user/gent/447/vsc44746/lidar_results/tree1-10m-002-3scans_point.txt"



data <- fread(input_file)
colnames(data)[1:3] <- c("X", "Y", "Z")
data <- data[, .(X, Y, Z)]

# Parameter initialization
D <- 0.01
cl_dist <- 0.01
max_d <- 0.10

# Add necessary columns
data[, `:=`(ID = .I, iter = -1)]

# First layer
lay <- data[Z <= min(Z) + 0.1]
data[lay$ID, iter := 1]

# Prepare data for hierarchical calculation
dat2 <- data[-lay$ID]

# Iterative layer calculation
i <- 2

while (nrow(dat2) > 0) {
  
  # Calculate maximum distance from each point to the previous layer
  dat2[, dist := Rfast::rowMaxs(FNN::knnx.dist(data = lay[, 1:3], query = dat2[, 1:3], algorithm = "kd_tree", k = 1), value = TRUE)]
  
  # Select the points for the next layer
  lay <- dat2[dist <= D]
  
  # If no points are found, select the nearest point to the already classified points as the new starting point
  if (nrow(lay) == 0) {
    lay <- dat2[which.min(FNN::knnx.dist(data = data[iter != -1, 1:3], query = dat2[, 1:3], algorithm = "kd_tree", k = 1))]
  }
  
  # Update iter value for the current layer points
  data[lay$ID, iter := i]
  
  # Remove the current layer points from the data
  dat2 <- dat2[!ID %in% lay$ID]
  i <- i + 1
}

#######################################################
#- clustering non connected components in each layer -#
#######################################################

# Perform clustering calculation for each layer
for (i in sort(unique(data$iter))) {
  
  # Get the indices of points in the current layer
  in_iter <- which(data$iter == i)
  
  if (length(in_iter) >= 2) {
    
    # Calculate the distance matrix and perform hierarchical clustering
    dist_matrix <- stats::dist(data[in_iter, .(X, Y, Z)])
    hclust_result <- fastcluster::hclust(dist_matrix, method = "single")
    clusters <- stats::cutree(hclust_result, h = cl_dist)
    
    # Store clustering results in data
    data[in_iter, cluster := clusters]

    
  } else {
    
    # If only one point, set it as a single cluster
    data[in_iter, cluster := 1]
  }
}
########################
#- Build the skeleton -#
########################

# Calculate the center of each cluster
cl <- data[, .(X = mean(X), Y = mean(Y), Z = mean(Z), iter = iter), by = .(iter, cluster)]
cl <- cl[, 3:6]

# Initialization
cl[, ':='(ID = 1:.N, done = 0)]
cl2 <- copy(cl)
root <- cl2[Z == min(Z)]

# Mark the root node as processed
cl[root$ID, done := 1]
cl2 <- cl2[-root$ID]

# Initialize the skeleton data table
skel <- data.table::data.table(matrix(ncol = 6))
setnames(skel, c("startX", "startY", "startZ", "endX", "endY", "endZ"))

# Build the skeleton
while (nrow(cl2) > 0) {
  
  # Record the coordinates of the current root node
  start <- root[, 1:3]
  
  # Calculate the distance of all unprocessed nodes to the current root node
  cl2[, dist := sqrt((X - root$X)^2 + (Y - root$Y)^2 + (Z - root$Z)^2)]
  
  # Find the next root node that meets the condition
  root <- cl2[dist == min(dist) & dist <= max_d & iter > root$iter]
  
  if (nrow(root) > 0) {
    # If a suitable root node is found, update it
    cl[root$ID, done := 1]
    cl2 <- cl2[!ID %in% root$ID]
    skel <- rbindlist(list(skel, data.table(start[, 1:3], root[, 1:3])), use.names = FALSE)
  } else {
    # If not found, select the lowest unprocessed node and connect it to the nearest processed node
    root <- cl2[which.min(iter)]
    root <- root[which.min(FNN::knnx.dist(data = cl[done == 1, 1:3], query = root[, 1:3], algorithm = "kd_tree", k = 1))]
    d <- sqrt((cl$X[cl$done == 1] - root$X)^2 + (cl$Y[cl$done == 1] - root$Y)^2 + (cl$Z[cl$done == 1] - root$Z)^2)
    done <- cl[done == 1]
    start <- done[d == min(d)]
    
    cl[root$ID, done := 1]
    cl2 <- cl2[ID != root$ID]
    skel <- rbindlist(list(skel, data.table(start[, 1:3], root[, 1:3])), use.names = FALSE)
  }
}

# Remove NA rows and calculate segment lengths
skel <- skel[-1,]
skel[, length := sqrt((startX - endX)^2 + (startY - endY)^2 + (startZ - endZ)^2)]
skel <- skel[length > 0]

######################
#- compute topology -#
######################
# Create segment ID
skel[, seg_ID := .I]  # Use .I for more efficient row number generation as seg_ID

# Assign bearer ID
skel[, bearer_ID := FNN::knnx.index(data = skel[, 4:6], query = skel[, 1:3], algorithm = "kd_tree", k = 1)]
skel[seg_ID == bearer_ID, bearer_ID := 0]  # If seg_ID is the same as bearer_ID, set bearer_ID to 0

# Calculate the total length borne by each segment
skel[, bear_length := length]  

# Iterate in descending order of seg_ID to calculate the total bear_length
for(s in skel[, rev(order(seg_ID))]) {
  # Find child nodes and sum bear_length
  skel[seg_ID == s, bear_length := bear_length + skel[bearer_ID == s, sum(bear_length)]]
}

# Identify axes）
skel[, axis_ID := 0]

# Select the starting segment as the initial axis
cur_seg <- skel[bearer_ID == 0]  # Starting segment (typically the base of the trunk)
cur_ID <- 1  # Current axis ID
cur_sec <- 1  # Current section

skel[bearer_ID == 0, axis_ID := cur_ID]  # Assign axis_ID to the starting segment


queue <- integer(0)  
while (any(skel$axis_ID == 0)) {
  
  skel[seg_ID == cur_seg$seg_ID, section := cur_sec]
  childs <- skel[bearer_ID == cur_seg$seg_ID]
  
  if (nrow(childs) == 1) {  # Only one child node -> belongs to the same axis
    skel[seg_ID == childs$seg_ID, axis_ID := cur_ID]
    cur_seg <- childs
  } else if (nrow(childs) > 1) {  # Multiple child nodes -> select the node bearing the longest structure
    max_child_idx <- which.max(childs$bear_length)
    skel[seg_ID == childs$seg_ID[max_child_idx], axis_ID := cur_ID]
    cur_seg <- childs[max_child_idx, ]
    
    queue <- c(queue, childs$seg_ID[-max_child_idx])
    cur_sec <- cur_sec + 1
  } else {  # No child nodes -> select the next node from the queue and increment cur_ID
    cur_ID <- cur_ID + 1
    cur_ID <- cur_ID + 1
    cur_seg <- skel[seg_ID == queue[1]]
    skel[seg_ID == cur_seg$seg_ID, axis_ID := cur_ID]
    
    queue <- queue[-1]
    
    cur_sec <- cur_sec + 1
  }
}


# ADD BRANCHING ORDER
cur_BO = skel[skel$axis_ID == 1] # axes of branching order 1
skel[axis_ID == 1,branching_order := 1]
BO = 2 # first branching order to detect
while(nrow(cur_BO)>0){
  # find all child axes of the axes of the curent BO
  child_axes=skel[bearer_ID %in% cur_BO$seg_ID &
                    !(axis_ID %in% unique(cur_BO$axis_ID)),c(axis_ID)]
  # add the new BO to the child axes
  skel[axis_ID %in% child_axes,branching_order := BO]
  # select the child axes for the next round
  cur_BO = skel[skel$axis_ID %in% child_axes]
  BO = BO+1
}

# Create the new_skeleton data table
new_skeleton <- skel[, .(
  startX, startY, startZ, endX, endY, endZ, 
  cyl_ID = seg_ID, parent_ID = bearer_ID, 
  extension_ID = 0, length, axis_ID, 
  segment_ID = section, node_ID = 0, branching_order
)]

# Optimize the segment_ID
new_skeleton[, segment_ID := max(cyl_ID), by = segment_ID]

# Preserve the original node_ID logic
new_skeleton[, node_ID := new_skeleton$segment_ID[which(
  new_skeleton$cyl_ID %in% parent_ID & new_skeleton$segment_ID != segment_ID
)], by = segment_ID]

# Calculate and refine extension_ID
new_skeleton[, extension_ID := FNN::knnx.index(data = new_skeleton[, .(startX, startY, startZ)], 
                                               query = new_skeleton[, .(endX, endY, endZ)], 
                                               algorithm = "kd_tree", k = 1)]

new_skeleton[cyl_ID == extension_ID | !cyl_ID %in% parent_ID, extension_ID := 0]


######################
#- Compute branch angles -#
######################

# Extract unique node_IDs
node_ID <- data.frame(new_skeleton$node_ID)
node<- unique(node_ID)

# Initialize result data frames
start_result <- data.frame() 
angle_result <- data.frame() 

# Get the number of unique nodes
number <- length(node$new_skeleton.node_ID)

# Iterate through each node starting from the second one
for (i in 2:number) {
  
  # Subset the new_skeleton data for the current node
  branch <- subset(new_skeleton, node_ID==node[i,]) 
  branch1 <- branch
  
  # Split the data into groups based on axis_ID
  grouped_data <- branch1 %>%
    group_by(axis_ID) %>%
    group_split()
  
  # If there are more than one group (i.e., branching occurs)
  if (length(grouped_data)>1){
    for (j in 2:length(grouped_data)) {
      
      # Extract the first group and the current group
      group3 <- grouped_data[[1]]     
      group4 <- grouped_data[[j]] 
      
      # Extract the first row from each group
      a <- group3[1,] 
      b <- group4[1,]
      
      # Compute the distance between the first points of both groups
      distance <- c(b-a) 
      
      # Adjust the current group by subtracting the distance
      group4 <- group4[-1, ]-distance
      
      # Calculate vectors AB and AC for angle computation
      A = unlist(group3[1,1:3])
      B = unlist(group3[nrow(group3), 4:6])
      C = unlist(group4[nrow(group4), 4:6])
      AB <- B - A
      AC <- C - A
      
      # Compute the dot product of vectors AB and AC
      dot_product <- sum(AB * AC)
      
      # # Compute vectors AB and AC
      norm_AB <- sqrt(sum(AB^2))
      norm_AC <- sqrt(sum(AC^2))
      
      # Convert the angle from radians to degrees (optional)
      angle_rad <- acos(dot_product / (norm_AB * norm_AC))
      
      # Append the results to the result data frames
      angle <- angle_rad * (180 / pi)
      
      start_result <- rbind(start_result, a)
      angle_result <- rbind(angle_result, angle)
      
      # Combine start_result and angle_result into the final result
      result <- cbind(start_result, angle_result)
    }
  }
}
result <- na.omit(result) 
result <- round(result, digits = 2)

# save automated extraction branch angles
write.table(result, file = output_file, row.names = FALSE, col.names = TRUE)
 2)

# save automated extraction branch angles
write.table(result, file = output_file, row.names = FALSE, col.names = TRUE)
