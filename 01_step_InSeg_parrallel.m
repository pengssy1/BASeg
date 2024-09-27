clear;clc;close all;
% Add path and import point cloud
addpath(genpath('A:\TreeQSM-master-241\src')); % Add TreeQSM path
addpath(genpath('C:\Users\xipeng\Downloads\skeletonization-master\cloudcontr_2_0')); %Add laplacian path
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
branch_size = accumarray(data2(:,1), data2(:,2), [], @max);
branch_number = unique(data1(:,4));

branch = [branch_number, branch_size];


data1(:, 5) = NaN;


for i = 1:size(branch, 1)
    order = branch(i, 1);
    size_value = branch(i, 2);
    
    
    data1(data1(:, 4) == order, 5) = size_value;
end


% Update branch segmentation data
[idx, ~] = knnsearch(data1(:, 1:3), data0(:, 1:3));

data0(:, 4) = data1(idx, 4);
data0(:, 5) = data1(idx, 5);

A = [data1; data0];

options.USING_POINT_RING = GS.USING_POINT_RING;  

% Get all unique order values
branch_orders = unique(A(:, 4));

% Initialize the matrix to store the results
final_result = [];

% computing to extract skeleton data
% for i = 1:length(branch_orders)
for i = 1:length(branch_orders)
    % current order 
    order = branch_orders(i);
    
    B = A(A(:, 4) == order, :);
    
    % branch size
    max_size = max(B(:, 5));
    
    % Set the downsampling resolution
    if max_size > 0.10
        grid_step = 0.04;
    elseif max_size > 0.02 && max_size <= 0.10
        grid_step = 0.02;
    else
        grid_step = 0.005;
    end
    
    % downsample each branch point cloud data
    ptCloud = pointCloud(B(:, 1:3));
    ptCloud_ds = pcdownsample(ptCloud, 'gridAverage', grid_step);
    
    C = ptCloud_ds.Location;

    D.pts = C;  
    D.npts = size(D.pts, 1);  
    [D.bbox, D.diameter] = GS.compute_bbox(D.pts);
    % Build local 1-ring       
    D.k_knn = GS.compute_k_knn(D.npts);
    D.rings = compute_point_point_ring(D.pts, D.k_knn, []);
    % Use the Laplace method to contract the point cloud
    [D.cpts, t, initWL, WC, sl] = contraction_by_mesh_laplacian(D, options);
    % save final skeleton data
    final_result = [final_result; D.cpts];
    
end

% save skeleton data
save('C:\Users\xipeng\Downloads\skeletonization-master\cloudcontr_2_0\result\tree_32_1.txt', 'final_result', '-ascii');







