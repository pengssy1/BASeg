clear;clc;close all;
% Add path and import point cloud
addpath(genpath('A:\TreeQSM-master-241\src')); % Add TreeQSM path
addpath(genpath('C:\Users\xipeng\Downloads\skeletonization-master\cloudcontr_2_0')); %A dd laplacian path
P = readmatrix("C:\Users\xipeng\Downloads\L1-Tree-main\Data\tree_A_032.txt"); % Read point cloud data in txt format

% Get branch segmentation input parameters
tic
inputs = define_input(P,1,1,1);

% Get initial variables and parameters
PatchDiam1 = inputs.PatchDiam1;
PatchDiam2Min = inputs.PatchDiam2Min;
PatchDiam2Max = inputs.PatchDiam2Max;
BallRad1 = inputs.BallRad1; 
BallRad2 = inputs.BallRad2; 

% Ensure that the point cloud data is 3D and double precision
P = double(P(:, 1:3));

% Branch point cloud segmentation
for h = 1:length(PatchDiam1)
  
  Inputs = inputs;
  Inputs.PatchDiam1 = PatchDiam1(h);
  Inputs.BallRad1 = BallRad1(h);
  
  % Generate initial cover sets
  cover1 = cover_sets(P, Inputs);
  
  % Determine tree sets and update neighbors
  [cover1, Base, Forb] = tree_sets(P, cover1, Inputs);
  
  % Generate initial segments
  segment1 = segments(cover1, Base, Forb);
  
  % Corrected segments
  segment1 = correct_segments(P, cover1, segment1, Inputs, 0, 1, 1);
  
  for i = 1:length(PatchDiam2Max)
    % Modify input parameters
    Inputs.PatchDiam2Max = PatchDiam2Max(i);
    Inputs.BallRad2 = BallRad2(i);

    for j = 1:length(PatchDiam2Min)
      % Modify input parameters
      Inputs.PatchDiam2Min = PatchDiam2Min(j);
      
      % Generate new cover sets
      RS = relative_size(P, cover1, segment1);
      cover2 = cover_sets(P, Inputs, RS);
      
      % Determine new tree sets and update neighbors
      [cover2, Base, Forb] = tree_sets(P, cover2, Inputs, segment1);
      
      % Generate new segments
      segment2 = segments(cover2, Base, Forb);
      
      % Modify segments and generate cylinders
      segment2 = correct_segments(P, cover2, segment2, Inputs, 1, 1, 0);
      cylinder = cylinders(P, cover2, segment2, Inputs);
    end
  end
end

data1 = [P, double(segment2.SegmentOfPoint)];
data0 = data1(data1(:, 4) == 0, :);
data2 = [double(cylinder.branch), double(cylinder.radius)];

% Delete the rows where the fourth column in data1 is 0
data1(data1(:, 4) == 0, :) = [];

% Extract the size of each branch
[unique_orders, ~, idx] = unique(data2(:, 1), 'stable');
first_branch_sizes = [unique_orders, accumarray(idx, data2(:, 2), [], @min)];

% Update branch segmentation data
result = nan(size(data1, 1), size(data1, 2) + 1);
for i = 1:size(data1, 1)
    current_order = data1(i, 4);
    size_value = first_branch_sizes(first_branch_sizes(:, 1) == current_order, 2);
    
    if ~isempty(size_value)
        result(i, :) = [data1(i, :), size_value];
    end
end

[idx, ~] = knnsearch(result(:, 1:3), data0(:, 1:3));

data0(:, 4) = result(idx, 4);
data0(:, 5) = result(idx, 5);

A = [result; data0];

options.USING_POINT_RING = GS.USING_POINT_RING;  % 设置选项

% Get all unique order values
branch_orders = unique(A(:, 4));

% Initialize the matrix to store the results
final_result = [];

% Use parallel computing to extract skeleton data
% parfor k = 1:length(branch_orders)
parfor k = 1:2
    order = branch_orders(k);
    
    % Retrieve the data of the current branch order
    B = A(A(:, 4) == order, :);

    % Initialize temporary variables to store results
    temp_result = [];

    % Downsample the data and select different resolutions according to the value of size
    for size_category = {'large', 'medium', 'small'}
        if strcmp(size_category, 'large')
            size_condition = B(:, 5) > 10;
            grid_step = 0.04;
        elseif strcmp(size_category, 'medium')
            size_condition = B(:, 5) > 2 & B(:, 5) <= 10;
            grid_step = 0.02;
        else
            size_condition = B(:, 5) <= 2;
            grid_step = 0.003;
        end
        
        % If there is point cloud data that meets the conditions
        if any(size_condition)
            % Extract data that meets the conditions
            data_to_downsample = B(size_condition, :);
            
            % Generate point cloud object and perform downsampling
            ptCloud = pointCloud(data_to_downsample(:, 1:3));
            ptCloud_ds = pcdownsample(ptCloud, 'gridAverage', grid_step);
            
            % Update the downsampled point cloud data
            downsampled_data = [ptCloud_ds.Location, data_to_downsample(1:size(ptCloud_ds.Location, 1), 4:5)];
            branch_data = downsampled_data(:, 1:3);

            % Process point cloud data
            C = struct();
            C.pts = branch_data;
            C.npts = size(C.pts, 1);
            [C.bbox, C.diameter] = GS.compute_bbox(C.pts);

            % Build local 1-ring
            C.k_knn = GS.compute_k_knn(C.npts);
            C.rings = compute_point_point_ring(C.pts, C.k_knn, []);

            % Use the Laplace method to contract the point cloud
            [C.cpts, ~, ~, ~, ~] = contraction_by_mesh_laplacian(C, options);

            % Add the current result to a temporary variable
            temp_result = [temp_result; C.cpts];
        end
    end
    
    % Add the temporary result to the final result
    final_result = [final_result; temp_result];
end
toc
% save skeleton data
save('C:\Users\xipeng\Downloads\skeletonization-master\cloudcontr_2_0\result\tree_32_test.txt', 'final_result', '-ascii');

