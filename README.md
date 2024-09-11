# InSeg

InSeg: A method to automatically extract branch angles from TLS data

How to implement it

1. Get tree skeleton data
Download TreeQSM 2.4.1  https://github.com/InverseTampere/TreeQSM
Download Laplacian algorithm  https://github.com/jjcao/skeletonization
Run 01_step_InSeg_parrallel.m in MATLAB, it will generate skeleton data of an individual tree. 
Replace the input file with an individual tree point cloud file in txt format in your local computer.

2. Automatically extract branch angles
2.1 Implement it in R
Run 02_angle_calculate in R, it will automatically extract branch angles from the skeleton data generated in the first step. 
The input path is the skeleton data in txt format generated in the previous step.

2.2 Implement it on HPC (batch submission script)
Upload 03_hpc to your HPC path

cd /data/gent/vo/000/gvo00074/Xipeng/creat # Enter the creat file on your HPC

python input_file_name.py # Extract the input tree point cloud file name, replace directory_path as your input file path on your HPC

python create_task.py # Generate corresponding run.r and job.sh files for all tree point clouds, replace run_path, driver_path, output_path as your file path on your HPC

cd /data/gent/vo/000/gvo00074/Xipeng/creat/run

bash run.sh # Batch submission script
