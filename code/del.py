import os
# the path of folder
folder_path = r'../code' 
file_list = os.listdir(folder_path)
# loop for delete .log and .trn
for file_name in file_list:
    if file_name.endswith('.log') or file_name.endswith('.trn'):
        file_path = os.path.join(folder_path, file_name)  
        os.remove(file_path) 