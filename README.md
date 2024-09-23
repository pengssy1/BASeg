# BASeg

# Retrieval of branch angles from terrestrial laser scanning data using an Instance Segmentation-based Contraction method

# How to implement BASeg

# 1. Download dependency packages

   TreeQSM 2.4.1  https://github.com/InverseTampere/TreeQSM
   
   Laplacian algorithm  https://github.com/jjcao/skeletonization

# 2. Get tree skeleton data

   Run 01_step_InSeg_parrallel.m in MATLAB.

   It will generate skeleton data of an individual tree. Replace the input file with an individual tree point cloud file in txt format in your local computer.

# 3. Automatically extract branch angles
   
   # 3.1 Implement it in R

   Run 02_angle_calculate in R.

   It will automatically extract branch angles from the skeleton data generated in the first step. The input path is the skeleton data in txt format generated in the previous step.

   # 3.2 Implement it on HPC (batch submission script)

   (1) Upload 03_hpc to your HPC path

   (2) Enter the creat file on your HPC
   
      cd /data/gent/vo/000/gvo00074/Xipeng/creat 

   (3) Extract the input tree point cloud file name, replace directory_path as your input file path on your HPC                 
   
      python input_file_name.py  

   (4) Generate corresponding run.r and job.sh files for all tree point clouds, replace run_path, driver_path, output_path as your file path on your HPC                   
   
      python create_task.py    

   (5) Enter run file on your HPC
   
      cd /data/gent/vo/000/gvo00074/Xipeng/creat/run

   (6) Batch submission script                               
   
      bash run.sh   
