import os

def read_txt_file(filename):
    with open(filename, 'r') as file:
        content = file.read()
    return content

def read_txt_to_list(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        lines = [line.strip() for line in lines]
    return lines

def append_line(filename, line):
    with open(filename, 'a') as file:
        file.write(line + '\n')

def copy_file(source, destination):
    with open(source, 'rb') as f_source:
        with open(destination, 'wb') as f_destination:
            # Read and write the file in chunks to conserve memory
            for chunk in iter(lambda: f_source.read(4096), b''):
                f_destination.write(chunk)
          
                
def create_task(site_id='line',
            run_path =r'C:/Users/xipeng/OneDrive - UGent/Desktop/create/run/',
            driver_path=r'./input_path/',output_path = r'./Outputs/'):
    
    task_dir = run_path + site_id +'/'
    
    if not os.path.exists(task_dir):
        os.mkdir(task_dir)
    
        copy_file('./job_template.sh',task_dir+'./job.sh')
        
        with open(task_dir+'./job.sh', 'r+') as file:
            lines = file.readlines()
            
            lines.insert(1, '#PBS -N ' + site_id + '\n')
            file.seek(0)
            
            file.writelines(lines)

        append_line(task_dir+'job.sh', 'Rscript ' + task_dir +"run.r")


        
        copy_file('./run.r',task_dir+'run.r')

        with open(task_dir+'run.r', 'r+') as file:
            lines = file.readlines()
            
            lines.insert(1, 'input_file =' + '"' + driver_path+site_id + '"\n')
            lines.insert(2, 'output_file=' + '"' + output_path + site_id + '.txt'+ '"\n')

            file.seek(0)
            
            file.writelines(lines)


if __name__ == '__main__':
        
    if not os.path.exists('../run.sh'):
        copy_file('./run.sh','./run/run.sh')

    site_info_file = 'site.txt'

    file_lines = read_txt_to_list(site_info_file)
    for line in file_lines:
        print(line)
        create_task(
            site_id=line,
            run_path =r'/data/gent/vo/000/gvo00074/Xipeng/creat/run/',
            driver_path=r'/data/gent/vo/000/gvo00074/Xipeng/creat/Inputs/',
            output_path=r'/data/gent/vo/000/gvo00074/Xipeng/creat/Outputs/')