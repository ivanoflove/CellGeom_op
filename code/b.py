import pandas as pd

# 读取CSV文件
data = pd.read_csv('../data/solution.csv', delimiter=',')



# 提取name和value列
name_values = data[['name', 'value']]

# 将name和value组合成name = value的形式
name_values['name_value'] = name_values['name'] + ' = ' + name_values['value'].astype(str)

# 选择name_value列并保存为txt文件
name_values['name_value'].to_csv('output.txt', index=False)