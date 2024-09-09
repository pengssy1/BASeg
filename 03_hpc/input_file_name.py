import os
directory_path = '/data/gent/vo/000/gvo00074/Xipeng/creat/Inputs'

filenames = os.listdir(directory_path)

output_file = 'site.txt'
with open(output_file, 'w') as file:
    for filename in filenames:
        file.write(filename + '\n')

print(f"Filenames have been written to {output_file}")
