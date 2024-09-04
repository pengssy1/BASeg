% InSet is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% ------------------------------------------------------------------------------
% INPUTS:
% fileName    Filename for point cloud data
%
% OUTPUTS:
% branch_angle.txt    The output file of the InSeg method, which records branch angles size and coordinates
% ------------------------------------------------------------------------------

clear;clc;close all;
% 1. import point cloud data
addpath(genpath('A:\TreeQSM-master-241\src'));
addpath(genpath('C:\Users\xipeng\Downloads\skeletonization-master\cloudcontr_2_0'));
P = readmatrix("C:\Users\xipeng\Downloads\L1-Tree-main\Data\tree_12.txt");

% 2. 获取输入参数
inputs = define_input(P,1,1,1);

% 3. 获取初始变量和参数
PatchDiam1 = inputs.PatchDiam1;
PatchDiam2Min = inputs.PatchDiam2Min;
PatchDiam2Max = inputs.PatchDiam2Max;
BallRad1 = inputs.BallRad1; 
BallRad2 = inputs.BallRad2; 

% 确保点云数据是3维且为双精度
P = double(P(:, 1:3));

% 4. 重建QSM
for h = 1:length(PatchDiam1)
  % 修改输入参数
  Inputs = inputs;
  Inputs.PatchDiam1 = PatchDiam1(h);
  Inputs.BallRad1 = BallRad1(h);
  
  % 生成初始的 cover sets
  cover1 = cover_sets(P, Inputs);
  
  % 确定 tree sets 并更新邻居
  [cover1, Base, Forb] = tree_sets(P, cover1, Inputs);
  
  % 生成初始 segments
  segment1 = segments(cover1, Base, Forb);
  
  % 修正 segments
  segment1 = correct_segments(P, cover1, segment1, Inputs, 0, 1, 1);
  
  for i = 1:length(PatchDiam2Max)
    % 修改输入参数
    Inputs.PatchDiam2Max = PatchDiam2Max(i);
    Inputs.BallRad2 = BallRad2(i);

    for j = 1:length(PatchDiam2Min)
      % 修改输入参数
      Inputs.PatchDiam2Min = PatchDiam2Min(j);
      
      % 生成新的 cover sets
      RS = relative_size(P, cover1, segment1);
      cover2 = cover_sets(P, Inputs, RS);
      
      % 确定新的 tree sets 并更新邻居
      [cover2, Base, Forb] = tree_sets(P, cover2, Inputs, segment1);
      
      % 生成新的 segments
      segment2 = segments(cover2, Base, Forb);
      
      % 修正 segments 并生成 cylinders
      segment2 = correct_segments(P, cover2, segment2, Inputs, 1, 1, 0);
      cylinder = cylinders(P, cover2, segment2, Inputs);
    end
  end
end

% 5. 处理 segment2.SegmentOfPoint 数据
data1 = [P, double(segment2.SegmentOfPoint)];
data2 = [double(cylinder.branch), double(cylinder.radius)];

% 删除 data1 中第四列为0的行
data1(data1(:, 4) == 0, :) = [];

% 提取每个 branch_order 的第一个 branch_size
[unique_orders, ~, idx] = unique(data2(:, 1), 'stable');
first_branch_sizes = [unique_orders, accumarray(idx, data2(:, 2), [], @min)];

% 6. 合并数据
result = nan(size(data1, 1), size(data1, 2) + 1);
for i = 1:size(data1, 1)
    current_order = data1(i, 4);
    size_value = first_branch_sizes(first_branch_sizes(:, 1) == current_order, 2);
    
    if ~isempty(size_value)
        result(i, :) = [data1(i, :), size_value];
    end
end

% 7. 保存结果
% save('C:\Users\xipeng\Downloads\skeletonization-master\cloudcontr_2_0\result.txt', 'result', '-ascii');


A = result;

options.USING_POINT_RING = GS.USING_POINT_RING;  % 设置选项

% 获取所有唯一的 order 值
unique_orders = unique(A(:, 4));

% 初始化存储结果的矩阵
final_result = [];

% 循环遍历每一个唯一的 order 值
for order = unique_orders'
    % 检索出当前 order 的数据
    B = A(A(:, 4) == order, :);

    % 对数据进行下采样，根据 size 的值选择不同的分辨率
    for size_category = {'large', 'medium', 'small'}
        if strcmp(size_category, 'large')
            size_condition = B(:, 5) > 10;
            grid_step = 0.04;
        elseif strcmp(size_category, 'medium')
            size_condition = B(:, 5) > 2 & B(:, 5) <= 10;
            grid_step = 0.02;
        else
            size_condition = B(:, 5) <= 2;
            grid_step = 0.005;
        end
        
        % 如果存在满足条件的点云数据
        if any(size_condition)
            % 提取满足条件的数据
            data_to_downsample = B(size_condition, :);
            
            % 生成点云对象并进行下采样
            ptCloud = pointCloud(data_to_downsample(:, 1:3));
            ptCloud_ds = pcdownsample(ptCloud, 'gridAverage', grid_step);
            
            % 更新下采样后的点云数据
            downsampled_data = [ptCloud_ds.Location, data_to_downsample(1:size(ptCloud_ds.Location, 1), 4:5)];
            branch_data = downsampled_data(:, 1:3);
            % Step 1: 读取点云数据
            tic
            C.pts = branch_data;  % 使用 readmatrix 读取 TXT 文件
            C.npts = size(C.pts, 1);  % 点的数量
            [C.bbox, C.diameter] = GS.compute_bbox(C.pts);
            toc

            % Step 2: 构建局部 1-ring
            tic
            C.k_knn = GS.compute_k_knn(C.npts);
            C.rings = compute_point_point_ring(C.pts, C.k_knn, []);
            toc

            % Step 3: 使用拉普拉斯方法收缩点云
            tic
            [C.cpts, t, initWL, WC, sl] = contraction_by_mesh_laplacian(C, options);
            toc

            % 将每个 order 的 result 数据保存到最终结果矩阵中
            final_result = [final_result; C.cpts];
        end
    end
end


% 保存最终结果到文件中
save('C:\Users\xipeng\Downloads\skeletonization-master\cloudcontr_2_0\skeleton_12.txt', 'final_result', '-ascii');




