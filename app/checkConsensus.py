# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 23:34:50 2020

@author: hill103
"""

"""
评估RNAhybrid和IntaRNA预测结果的一致性
一致性判断基于checkDuplicate函数（确认IntaRNA结果中的重复entries）改写而成
经过考虑，offset bases threshold设为2，即前后可以相差2个bases
"""


import pandas as pd
from time import time
from tqdm import tqdm
import sys, os
import collections
from getopt import getopt

start_time = time()


###############Consensus函数#####################################################
def checkConsensus(dataframe, offset):
    '''检查所有记录，返回重复的、一致的entries的index
    IntaRNA的重复entries已经delete，因此不必考虑重复entries的影响
    也因此consensus记录必定是一一对应，不会出现一对多现象
    '''
    
    def checkEntries(dataframe, offset):
        '''检查某个tRF和某个transcript的所有记录，是否属于一致
        一致的定义为：target起止位置在+/-offset个base之内
        dataframe只包含1个tRF和1个transcript
        返回一致entries的index
        '''
    
        def checkSite(start1, start2, end1, end2, offset):
            '''输入2个interaction site的开始和结束位置
            判断是否属于一致
            '''
        
            if abs(start1-start2)<=offset and abs(end1-end2)<=offset:
                return True
            else:
                return False
    
        if dataframe.shape[0] == 1:
            return None
    
        output = []
        # 从第2条记录开始
        for i in range(1, dataframe.shape[0], 1):
            # 确认与排在它前面的所有记录，是否存在一致
            for j in range(i):
                if checkSite(dataframe.at[dataframe.index[j], 'Start_Target'],
                             dataframe.at[dataframe.index[i], 'Start_Target'],
                             dataframe.at[dataframe.index[j], 'End_Target'],
                             dataframe.at[dataframe.index[i], 'End_Target'], offset):
                    # 配对二者的index都加入list
                    output.append(dataframe.index[i])
                    output.append(dataframe.index[j])
                    # 只要配对上了，该记录就可以保留，无需考虑它最终能配对上多少个
                    break
        
        # 删除output中的重复项，确认其数目是否为偶数
        # 注意保留原始list，因为set会自动排序
        # 原始list中是相邻两两配对的，利于网页展示
        if len(set(output)) % 2 == 1:
            # 存在一对多情况
            find_ids = [item for item, count in collections.Counter(output).items() if count > 1]
            print('########################################')
            print('Special cases exist for tRF "{}" and transcript "{}"'.format(
                    dataframe.iloc[0]['tRF_ID'], dataframe.iloc[0]['Transcript_ID']))
            for idx in find_ids:
                # 输出特例
                print(dataframe.at[idx, 'Demo'])
            print('########################################')
            
        return output
    
    
    # 先对dataframe进行排序
    dataframe.sort_values(['Transcript_ID'], inplace=True)
    
    to_del = []
    
    this_enst = dataframe.at[dataframe.index[0], 'Transcript_ID']
    ind_for_check = []
    
    # at/loc都是按index索引row，iloc按位置
    for i in tqdm(dataframe.index):        
        if dataframe.at[i, 'Transcript_ID'] == this_enst:
            # 同一个transcript的记录，加入比较
            ind_for_check.append(i)
        else:
            # 新的transcript的开始，之前的可以进行比较了
            tmp_result = checkEntries(dataframe.loc[ind_for_check], offset)
            if not tmp_result is None:
                if len(tmp_result) > 0:
                    '''
                    print('Total {:d} dumplicated entries for tRF "{}" and transcript "{}"'
                          .format(len(tmp_result), dataframe.at[tmp_result[0], 'tRF_ID'], this_enst))
                    '''
                    to_del += tmp_result
                    
            this_enst = dataframe.at[i, 'Transcript_ID']
            ind_for_check.clear()
            ind_for_check.append(i)
    
    # 循环结束时，最后一组序列需要检查
    if len(ind_for_check) > 0:
        tmp_result = checkEntries(dataframe.loc[ind_for_check], offset)
        if not tmp_result is None:
            if len(tmp_result) > 0:
                '''
                print('Total {:d} dumplicated entries for tRF "{}" and transcript "{}"'
                      .format(len(tmp_result), dataframe.at[tmp_result[0], 'tRF_ID'], this_enst))
                '''
                to_del += tmp_result
    
    # 最后返回to_del索引中的rows
    return dataframe.loc[to_del]
    

###############主函数#####################################################
def usage():
    '''对主函数进行简介
    '''
    print('''
python checkConsensus.py [option][value]...
    -h or --help        print this help messages
    -r or --rnahybrid   CSV files of parsed RNAhybrid results
    -i or --intarna     CSV files of parsed IntaRNA results
    -o or --outputpath  absolute or relative path for CSV file of consensus results. If ignored, the current path will be used
''')

# 如果没有任何参数，显示提示信息，并退出
if len(sys.argv) == 1:
    print('-h or --help for detail')
    sys.exit(1)


# 定义命令行参数
# 短选项名后的冒号(:)表示该选项必须有附加的参数
# 长选项名后的等号(=)表示该选项必须有附加的参数。
shortargs = 'hr:i:o:'
longargs = ['help', 'rnahybrid=', 'intarna=', 'outputpath=']

# 解析命令行参数
# sys.argv[0]为python脚本名，后续全为参数
# opts为分析出的参数信息，args为不符合格式信息的剩余参数
opts, args = getopt(sys.argv[1:], shortargs, longargs)


# 如果存在不符合格式信息的剩余参数，显示提示信息，并退出
if args:
    print('Invalid options exist!')
    print('-h or --help for detail')
    sys.exit(1)

   
# 定义dict类型的参数集，使得算法更稳健
paramdict = {'rnahybrid_file':None, 'intarna_file':None, 'output_path':os.getcwd()}
 
for opt,val in opts:
    if opt in ('-h', '--help'):
        usage()
        sys.exit(1)
        
    if opt in ('-r', '--rnahybrid'):
        if not os.path.isfile(val):
            # 输入不是一个确实存在的文件名
            raise Exception('Invalid input RNAhybrid CSV file!')
        # 采用realpath函数，获得真实绝对路径
        paramdict['rnahybrid_file'] = os.path.realpath(val)
        continue
    
    if opt in ('-i', '--intarna'):
        if not os.path.isfile(val):
            # 输入不是一个确实存在的文件名
            raise Exception('Invalid input IntaRNA CSV file!')
        # 采用realpath函数，获得真实绝对路径
        paramdict['intarna_file'] = os.path.realpath(val)
        continue
    
    if opt in ('-o', '--outputpath'):
        if not os.path.isdir(val):
            # 输入不是一个确实存在的文件夹
            raise Exception('Invalid or not existed directory!')
        # 采用realpath函数，获得真实绝对路径
        paramdict['output_path'] = os.path.realpath(val)
        continue

# 检查参数是否齐全
for k,v in paramdict.items():
    if v is None:
        raise Exception('Option "{}" is missing!'.format(k))

print('checking Consensus interactions options:')
for k,v in paramdict.items():
    print('{}: {}'.format(k, v))


# ---------------------Needed Files------------------------------
output_path = paramdict['output_path']
rnahybrid_file = paramdict['rnahybrid_file']
intarna_file = paramdict['intarna_file']
output_file = os.path.join(output_path, 'consensus_results.csv')


# Read all results at once, it will cost huge RAM

# Read RNAhybrid results
rnahybrid_df = pd.read_csv(rnahybrid_file)
print('Read {:,} RNAhybrid entries'.format(rnahybrid_df.shape[0]))

intarna_df = pd.read_csv(intarna_file)
print('Read {:,} IntaRNA entries'.format(intarna_df.shape[0]))

# 所有tRFs
trfs = set(rnahybrid_df.tRF_ID)
print('total {:d} tRFs'.format(len(trfs)))


# 统计结果的dataframe
status = pd.DataFrame(columns=['tRF_ID', '#RNAhybrid', '#IntaRNA',
                               '#Consensus', 'Pr_RNAhybrid', 'Pr_IntaRNA'])

# 遍历所有tRFs
for ind, trf in enumerate(trfs):
    
    print('Currently checking tRF "{}"'.format(trf))
    
    # Check Consensus
    # the columns of RNAhybrid and IntaRNA are already setted to the same order
    this_rnahybrid = rnahybrid_df.loc[rnahybrid_df.tRF_ID==trf]
    this_intarna = intarna_df.loc[intarna_df.tRF_ID==trf]
    
    output = checkConsensus(pd.concat([this_rnahybrid, this_intarna], ignore_index=True), 2)
    print('Finally get {:,} consensus entries'.format(output.shape[0]))
    print()
    
  
    # 数据存入CSV文件
    # if file does not exist write header
    if ind==0:
       output.to_csv(output_file, header=True, index=False)
    else: # else it exists so append without writing the header
       output.to_csv(output_file, mode='a', header=False, index=False)
    
    
    # less efficient
    status = status.append({'tRF_ID': trf,
                            '#RNAhybrid': this_rnahybrid.shape[0], 
                            '#IntaRNA': this_intarna.shape[0],
                            '#Consensus': output.shape[0],
                            'Pr_RNAhybrid': output.shape[0] / 2.0 / this_rnahybrid.shape[0],
                            'Pr_IntaRNA': output.shape[0] / 2.0 / this_intarna.shape[0]},
                           ignore_index=True)

# 保存统计结果
status.to_csv(os.path.join(output_path, 'tRF_level_consensus_stats.csv'), index=False)
    
print('All tRFs entries checked. Elapsed time: {:.2f} hours'.format(
        (time()-start_time)/3600.0))
