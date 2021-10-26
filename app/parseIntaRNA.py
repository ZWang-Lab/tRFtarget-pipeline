# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 02:59:10 2020

@author: hill103
"""

'''
Parse result from IntaRNA, and save into a MySQL database
本次自定义程序，生成interaction的示意图。经过与IntaRNA画出的示意图进行比对，本人程序生成的示意图均正确
函数getDemo，getMaxHitDP，checkDuplicate也用于解析RNAhybrid的结果
'''


import pandas as pd
from tqdm import tqdm
import re
from math import floor
from time import time
import sys, os
from getopt import getopt
import numpy as np


begin_time = time()


# ---------------------数据库相关-------------------------------
# 先保存成CSV文件


# ---------------------结果解析相关函数---------------------------
def parseTranID(tran_id):
    # 解析transcript的id，返回ENST编码和ENSG编码，不带版本号
    # transcript的id号为'|'分隔的多个字符串
    # 第1个为ENST编码，第2个为ENSG编码
    # 版本号为'.'分隔
    tmp_id = tran_id.strip().split('|')
    # return tmp_id[0].split('.')[0], tmp_id[1].split('.')[0]
    # 除human和mouse外，其余物种ENST和ENSG编码没有版本号
    return tmp_id[0], tmp_id[1]


def getArea(tran_id, start1):
    '''根据interaction在mRNA上的**起始位置**，以及mRNA的区域起止说明
    判断绑定点位于哪一个区域
    '''
    # 最后一个为空字符串，倒数2-4为'UTR5','CDS','UTR3'起止位置，但不一定3个部分都有
    tmp_id = tran_id.strip().split('|')
    tmp_part_dict = {}
    for i in range(-2, -5, -1):
        tmp_list = tmp_id[i].split(':')
        if len(tmp_list) == 2:
            tmp_part_dict[tmp_list[0]] = tmp_list[1]
    tmp_pos = int(start1)
    for k, v in tmp_part_dict.items():
        # 分隔起止位置
        tmp = v.split('-')
        if tmp_pos>=int(tmp[0]) and tmp_pos<=int(tmp[1]):
            return k


def getMaxHitLen(demo):
    '''给出一个匹配的示意图，找到连续匹配上的bases的最大长度
    G-U不稳定匹配也认为是匹配
    示意图共9行，其中第5行是匹配符号
    '''
    tmp = demo.split('\n')
    matched = tmp[4].strip().split()
    max_len = 0
    for seq in matched:
        if max_len < len(seq):
            max_len = len(seq)
    return max_len  


def getMaxHitDP(demo, max_len):
    '''给出一个匹配的示意图，找到最长的连续匹配上的序列
    G-U不稳定匹配也认为是匹配
    示意图共9行，其中第5行是匹配符号
    返回格式为target seq&query seq，并且二者都是5'->3'方向
    max_len为最长的连续匹配长度
    多个序列用'|'隔开
    '''
    tmp = demo.split('\n')
    sub_seq = [seq for seq in tmp[4].strip().split() if len(seq)==max_len]
    result = []

    for seq in sub_seq:
        # 确定序列在demo中的起始位置
        index = tmp[4].index(seq)
        # 用'&'连接target和query，并且query反向
        result.append(tmp[3][index:index+max_len]+'&'+tmp[5][index:index+max_len][::-1])
            
    return '|'.join(result)


def checkDuplicate(dataframe):
    '''检查所有记录，返回重复的、需要删除的entries的index
    dataframe已经经过排序
    '''
    
    def checkEntries(dataframe):
        '''检查某个tRF和某个transcript的所有记录，是否属于重复
        重复的定义为：target起止位置在+/-8个base之内
        dataframe只包含1个tRF和1个transcript，并且按MFE从低到高排序
        返回需要删除的重复entries的index
        重复entries只保留energy最低的记录
        '''
    
        def checkSite(start1, start2, end1, end2):
            '''输入2个interaction site的开始和结束位置
            判断是否属于重复
            '''
        
            if abs(start1-start2)<=8 and abs(end1-end2)<=5:
                return True
            elif abs(start1-start2)<=5 and abs(end1-end2)<=8:
                return True
            else:
                return False
    
        if dataframe.shape[0] == 1:
            return None
    
        output = []
        # 从第2条记录开始
        for i in range(1, dataframe.shape[0], 1):
            # 确认与排在它前面的所有记录，是否存在重复
            for j in range(i):
                if checkSite(dataframe.at[dataframe.index[j], 'Start_Target'],
                             dataframe.at[dataframe.index[i], 'Start_Target'],
                             dataframe.at[dataframe.index[j], 'End_Target'],
                             dataframe.at[dataframe.index[i], 'End_Target']):
                    output.append(dataframe.index[i])
                    break
            
        return output
    
    
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
            tmp_result = checkEntries(dataframe.loc[ind_for_check])
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
        tmp_result = checkEntries(dataframe.loc[ind_for_check])
        if not tmp_result is None:
            if len(tmp_result) > 0:
                '''
                print('Total {:d} dumplicated entries for tRF "{}" and transcript "{}"'
                      .format(len(tmp_result), dataframe.at[tmp_result[0], 'tRF_ID'], this_enst))
                '''
                to_del += tmp_result
            
    return to_del


# ---------------------Demo生成相关函数---------------------------
def getDemo(full_seq1, start1, end1, full_seq2, start2, end2, subseqDP, hybridDP):
    '''生成interaction的plain text示意图
    1是指target RNA
    2是指query tRF
    subseqDP是用'&'连接的interaction局部序列
    hybridDP是用'&'连接的interaction局部匹配模式图(dot-bracket notation)
    注意：1.start和end是以**1**为起点的，需要包含end，因此起止为[start1-1:end1]
    2.匹配示意图为target在上，tRF在下，并且tRF需要进行反向(3'->5')
    3.所有T变成U
    '''
    
    def getInteractDemo(seq1, seq2, match1, match2):
        '''构建interaction区域示意图
        tRF已经经过反转，从左至右处理即可
        '''
    
        def addNote(base1, base2):
            '''根据匹配情况，确定匹配符号
            '|'表示匹配，':'表示G-U不稳定匹配
            '''
            if (base1, base2) in set([('A','U'), ('U','A'), ('C','G'), ('G','C')]):
                return '|'
            elif (base1, base2) in set([('G','U'), ('U','G')]):
                return ':'
            else:
                raise Exception('Invalid complementary base pair "{}"-"{}"!'.format(base1, base2))
        
    
        # 5个空白行
        line1 = '' # target中非匹配base
        line2 = '' # target中匹配base
        line3 = '' # 匹配符号，'|'表示匹配，':'表示G-U不稳定匹配
        line4 = '' # tRF中匹配base
        line5 = '' # tRF中非匹配base
        # tRF和target分别使用index
        ind1 = 0
        ind2 = 0
    
        while (ind1<len(seq1)) and (ind2<len(seq2)):
            # 共4种符号pattern
            # 第1种，匹配
            if match1[ind1]=='(' and match2[ind2]==')':
                line1 += ' '
                line2 += seq1[ind1]
                line3 += addNote(seq1[ind1], seq2[ind2])
                line4 += seq2[ind2]
                line5 += ' '
                ind1 += 1
                ind2 += 1
            # 第2种，都不匹配
            elif match1[ind1]=='.' and match2[ind2]=='.':
                line1 += seq1[ind1]
                line2 += ' '
                line3 += ' '
                line4 += ' '
                line5 += seq2[ind2]
                ind1 += 1
                ind2 += 1
            # 第3种，target匹配，tRF不匹配
            # 此时target处出现空格，跳过1格
            elif match1[ind1]=='(' and match2[ind2]=='.':
                line1 += ' '
                line2 += ' '
                line3 += ' '
                line4 += ' '
                line5 += seq2[ind2]
                ind2 += 1
            # 第4种，target不匹配，tRF匹配
            # 此时tRF处出现空格，跳过1格
            elif match1[ind1]=='.' and match2[ind2]==')':
                line1 += seq1[ind1]
                line2 += ' '
                line3 += ' '
                line4 += ' '
                line5 += ' '
                ind1 += 1
            else:
                raise Exception('Exist unrecognized dot-bracket notation pattern!')
    
        # 确认示意图正确性
        assert(ind1==len(seq1))
        assert(ind2==len(seq2))
        
        return [line1, line2, line3, line4, line5]
    
    def getAddSeq(full_seq, ind, part, seq_type):
        '''生成interaction区域外部的序列示意图
        格式为3个base+3省略号+5个base
        part用于指示是5' head或3' tail区域
        seq_tytpe用于指示是target RNA或query tRF
        ind以**1**为起点
        输入tRF并未反转
        '''
        if part == 'head':
            # 5' head
            # 判断是否需要省略号
            if ind > 8:
                add_seq = full_seq[:3] + '...' + full_seq[ind-6:ind-1]
            else:
                add_seq = full_seq[:ind-1]
        else:
            # 3' tail区域
            # 首先判断是否需要省略号
            if len(full_seq)-ind > 8:
                add_seq = full_seq[ind:ind+5] + '...' + full_seq[-3:]
            else:
                add_seq = full_seq[ind:]
    
        # 加上5'或3'标记
        if seq_type == 'target':
            if part == 'head':
                return '5\'-'+add_seq
            else:
                return add_seq+'-3\''
        else:
            # 反转tRF
            if part == 'head':
                return add_seq[::-1]+'-5\''
            else:
                return '3\'-'+add_seq[::-1]
    
    def getNoteLine(start, end, len_pre, len_mid, len_suf, len_full, seq_type):
        '''确定line2或line8，包括指示序列位置的'|'和空格
        start和end以1为起点
        '''
        # 前缀
        if seq_type == 'target':
            # line2 for target, 5' head
            if start == 1:
                # interaction起始标记就是start
                line = ' '*3
            elif len_pre < 5:
                # interaction外部序列<2，为避免1和2挨在一起，无需起始指示
                line = ' '*len_pre
            else:
                line = ' '*3 + '|' + ' '*(len_pre-4)
        else:
            # line8 for query, 3' tail
            if end == len_full:
                # interaction终止就是序列end
                line = ' '*3
            elif len_pre < 6:
                # interaction外部序列<3，为避免数字重叠，无需首尾指示
                line = ' '*len_pre
            else:
                line = ' '*3 + '|' + ' '*(len_pre-4)
            
        # 加上interaction段
        line += ('|' + ' '*(len_mid-2) + '|')
    
        # 后缀
        if seq_type == 'target':
            # line2 for target, 3' tail
            # 如果interaction终止就是序列end，无需任何操作
            if end == len_full:
                line += ' '*3
            elif len_suf < 10:
                # interaction外部序列<7，为避免数字重叠，无需首尾指示
                line += ' '*len_suf
            else:
                line += (' '*(len_suf-4) + '|' + ' '*3)
        else:
            # line8 for query, 5' head
            # 如果interaction起点就是序列start，无需任何操作
            if start == 1:
                line += ' '*3
            elif len_suf < 5:
                # interaction外部序列<2，为避免1和2挨在一起，无需起始指示
                line += ' '*len_suf
            else:
                line += (' '*(len_suf-4) + '|' + ' '*3)
    
        return line

    def getNoteNum(line, start, end, len_full, seq_type):
        '''根据line2或line8，'|'的位置
        确定line1或line9，数字的位置
        暂不考虑数字之间是否会出现间隔太小导致覆盖的问题
        默认不会出现覆盖的现象
        '''
    
        def strReplace(string, index, insert_str):
            '''在string的index处，将同样长度的子字符串替换成insert_str
            字符串长度不变
            根据insert_str的长度进行自适应调整，使得insert_str在index处居中
            不考虑溢出，即string在index处留下的空间小于insert_str的长度
            '''
            tmp_len = len(insert_str)
            # index处之前需要被替换的字符数，即向前偏移量
            pre_len = floor(tmp_len/2)
            # index处会占用一个字符
            # index处之后，位置偏移量
            suf_len = tmp_len - pre_len
            return string[:index-pre_len] + insert_str + string[index+suf_len:]

        # 确定line中'|'的位置，注意使用反斜杠转义
        indices = [m.start() for m in re.finditer('\|', line)]
        # 确定对应的数字
        # 顺序为1, start, end, len_full
        nums = [1, start, end, len_full]
        if start == 1:
            del nums[0]
        if end == len_full:
            del nums[-1]
        # 如果此时nums仍然比竖线多，说明为避免数字重叠，序列首部或者末尾的指示被取消了
        if start == 2:
            del nums[0]
        
        if len(nums) > len(indices):
            del nums[-1]
    
        assert(len(nums)==len(indices))
    
        # query序列需要反向
        if seq_type == 'query':
            nums.reverse()
    
        # 采用字符替换的方法确定输出
        new_line = line
        for ind, num in zip(indices, nums):
            new_line = strReplace(new_line, ind, str(num))
        return new_line

    # 按'&'分隔target和tRF
    sub_seq1, sub_seq2 = subseqDP.split('&')
    sub_match1, sub_match2 = hybridDP.split('&')
    # T变成U
    tmp_full1 = full_seq1.replace('T', 'U')
    sub_seq1 = sub_seq1.replace('T', 'U')
    tmp_full2 = full_seq2.replace('T', 'U')
    sub_seq2 = sub_seq2.replace('T', 'U')
    
    # 构建interaction区域示意图，输入反转后的tRF
    # 输出为line3-line7的中间部分
    # interaction区域所有line长度一致
    demo = getInteractDemo(sub_seq1, sub_seq2[::-1], sub_match1, sub_match2[::-1])
    
    # interaction前后的序列示意图
    line3_pre = getAddSeq(tmp_full1, start1, 'head', 'target')
    line3_suf = getAddSeq(tmp_full1, end1, 'tail', 'target')    
    line7_pre = getAddSeq(tmp_full2, end2, 'tail', 'query')
    line7_suf = getAddSeq(tmp_full2, start2, 'head', 'query')
    
    # 确定line2，即target位置标记的竖线位置
    line2 = getNoteLine(start1, end1, len(line3_pre), len(demo[0]),
                        len(line3_suf), len(tmp_full1), 'target')
    # 确定line1，即target的位置数字
    line1 = getNoteNum(line2, start1, end1, len(tmp_full1), 'target')
    
    # 确定line8，即query位置标记的竖线位置
    line8 = getNoteLine(start2, end2, len(line7_pre), len(demo[0]),
                        len(line7_suf), len(tmp_full2), 'query')
    # 确定line9，即query的位置数字
    line9 = getNoteNum(line8, start2, end2, len(tmp_full2), 'query')
    
    # 前缀序列补齐，后缀不齐不影响
    # 看line3_pre和line7_pre谁长
    diff = abs(len(line3_pre) - len(line7_pre))
    if len(line3_pre) > len(line7_pre):
        # line7,8,9前面需要补空格
        line7_pre = ' '*(diff) + line7_pre
        line8 = ' '*(diff) + line8
        line9 = ' '*(diff) + line9
    elif len(line3_pre) < len(line7_pre):
        # line1,2,3前面需要补空格
        line1 = ' '*(diff) + line1
        line2 = ' '*(diff) + line2
        line3_pre = ' '*(diff) + line3_pre
        
    # 合并所有line
    tmp_len = len(line3_pre)
    # interaction区域为line3-line7的中间部分
    # 5'和3'标记可能加至不同的line上
    # line3/4，前缀
    if start1 == 1:
        demo[0] = ' '*tmp_len + demo[0]
        demo[1] = line3_pre + demo[1]
    else:
        demo[0] = line3_pre + demo[0]
        demo[1] = ' '*tmp_len + demo[1]
    # line3/4，后缀
    if end1 == len(tmp_full1):
        demo[1] = demo[1] + line3_suf
    else:
        demo[0] = demo[0] + line3_suf
    # line5
    demo[2] = ' '*tmp_len + demo[2]
    # line6/7，前缀
    if end2 == len(tmp_full2):
        demo[3] = line7_pre + demo[3]
        demo[4] = ' '*tmp_len + demo[4]
    else:
        demo[3] = ' '*tmp_len + demo[3]
        demo[4] = line7_pre + demo[4]
    # line6/7，后缀
    if start2 == 1:
        demo[3] = demo[3] + line7_suf
    else:
        demo[4] = demo[4] + line7_suf
        
    demo = [line1, line2] + demo + [line8, line9]
    
    # 最后在line4和line6前面加上target和tRF字样，占8个字符
    # 同时去除line1,2,8,9末尾存在的2个空格(对应的是3'和5'字符)
    for i in range(len(demo)):
        if i == 3:
            demo[i] = 'target  ' + demo[i].rstrip()
        elif i == 5:
            demo[i] = 'tRF     ' + demo[i].rstrip()
        else:
            demo[i] = ' '*8 + demo[i].rstrip()
    
    return '\n'.join(demo)



#############主函数#####################################################################
def usage():
    '''对主函数进行简介
    '''
    print('''
python parseIntaRNA.py [option][value]...
    -h or --help        print this help messages
    -d or --directory   directory of IntaRNA output file. The parsed results will also be saved in the same directory. support absolute or relative path
''')

# 如果没有任何参数，显示提示信息，并退出
if len(sys.argv) == 1:
    print('-h or --help for detail')
    sys.exit(1)


# 定义命令行参数
# 短选项名后的冒号(:)表示该选项必须有附加的参数
# 长选项名后的等号(=)表示该选项必须有附加的参数。
shortargs = 'hd:'
longargs = ['help', 'directory=']

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
paramdict = {'output_path':os.getcwd()}
 
for opt,val in opts:
    if opt in ('-h', '--help'):
        usage()
        sys.exit(1)
        
    if opt in ('-d', '--directory'):
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

print('IntaRNA results parsing options:')
for k,v in paramdict.items():
    print('{}: {}'.format(k, v))


# ---------------------Needed Files------------------------------
output_path = paramdict['output_path']
intarna_file = os.path.join(output_path, 'intarna_results.csv')
trf_info_file = os.path.join(output_path, 'trfs_info.csv')
tran_info_file = os.path.join(output_path, 'transcripts_info.csv')
output_file = os.path.join(output_path, 'parsed_intarna_results.csv')


# ---------------------准备附加信息------------------------------

# 1.tRF信息，用tRF ID索引
# 读取csv文件，并保存成dict格式，并且'T'变成'U
tmp_data = pd.read_csv(trf_info_file)
tRF_seq = {}
for i in tmp_data.index:
    tRF_seq[tmp_data.at[i, 'tRF_ID']] = tmp_data.at[i, 'tRF_Seq'].replace('T', 'U')
print('Total {:,} tRF sequences'.format(len(tRF_seq)))
del tmp_data


# 2.gene symbol信息dict gene_symbol，用gene ensembl id索引
# skipped
 

# 3.transcript type和name信息rna_infos，用transcript ensembl id索引
# skipped


# 4.transcript 序列信息rna_seq，用transcript ensembl id索引
# 读取csv文件，并保存成dict格式，并且'T'变成'U
tmp_data = pd.read_csv(tran_info_file)
rna_seq = {}
for i in tmp_data.index:
    rna_seq[tmp_data.at[i, 'Trans_ID']] = tmp_data.at[i, 'Trans_Seq'].replace('T', 'U')
print('Total {:,} transcripts'.format(len(rna_seq)))
del tmp_data

    
# ---------------------解析intaRNA结果------------------------------
# 读入CSV文件
start_time = time()
inta_result = pd.read_csv(intarna_file, sep=';', dtype={'id2':object})
print('All entries loaded. Elapsed time: {:.2f} minutes'.format(
        (time()-start_time)/60.0))

# 需要增加的新信息
# 注意：empty string不是null
inta_result['Max_Hit_Len'] = np.nan
inta_result['Demo'] = ''
inta_result['Max_Hit_DP'] = ''
inta_result['P_Val'] = np.nan
# 增加Tool
inta_result['Tool'] = 'IntaRNA'
        
# 列名重命名
inta_result.rename(columns={'hybridDP': 'HybridDP', 'subseqDP': 'SubseqDP',
                     'id2': 'tRF_ID', 'id1': 'Transcript_ID', 'E': 'MFE',
                     'start1': 'Start_Target', 'end1': 'End_Target',
                     'start2': 'Start_tRF', 'end2': 'End_tRF'}, inplace=True)

# inta_result.info(memory_usage='deep') # dataframe占用内存


# 确认重复entries
# 对所有entries进行排序
start_time = time()
inta_result.sort_values(['tRF_ID', 'Transcript_ID', 'MFE'],
                        ascending=[True, True, True], inplace=True)

# When index is unique, pandas use a hashtable to map key to value O(1)
# When index is non-unique and sorted, pandas use binary search O(logN)
# When index is non-unique and non-sorted, pandas need to check all the keys in the index O(N)
print('Sorting dataframe completed. Elapsed time: {:.2f} hours'.format((time()-start_time)/3600.0))

# 获得重复entries index的list
print('Begin checking duplicates')
start_time = time()
to_del = checkDuplicate(inta_result)
print('Total {:,} duplicates need to be deleted'.format(len(to_del)))

'''
# 删除重复entries by index
if len(to_del) > 0:
    inta_result.drop(to_del, inplace=True)
'''

# drop大量rows会非常慢
# 反过来从中抽取需要的rows
if len(to_del) > 0:
    need_iloc = []
    to_del_set = set(to_del)

    for i in inta_result.index:
        if i in to_del_set:
            need_iloc.append(False)
        else:
            need_iloc.append(True)

    inta_result2 = inta_result.iloc[need_iloc]
else:
    inta_result2 = inta_result

# Remain entries
print('Remain {:,} entries after delete duplicated entries'.format(inta_result2.shape[0]))
print('Elapsed time: {:.2f} hours'.format((time()-start_time)/3600.0))

del inta_result, to_del
if 'to_del_set' in locals():
    del to_del_set, need_iloc

# inta_result2.info(memory_usage='deep') # dataframe占用内存


# 解析获得其余的features
print('Start parsing the intaRNA result to get the rest features...')
for i in tqdm(inta_result2.index):
    
    inta_result2.at[i, 'Demo'] = getDemo(rna_seq[inta_result2.at[i, 'Transcript_ID'].strip()],
                  inta_result2.at[i, 'Start_Target'], inta_result2.at[i, 'End_Target'],
                  tRF_seq[inta_result2.at[i, 'tRF_ID'].strip()],
                  inta_result2.at[i, 'Start_tRF'], inta_result2.at[i, 'End_tRF'],
                  inta_result2.at[i, 'SubseqDP'], inta_result2.at[i, 'HybridDP'])
    
    inta_result2.at[i, 'Max_Hit_Len'] = getMaxHitLen(inta_result2.at[i, 'Demo'])
    inta_result2.at[i, 'Max_Hit_DP'] = getMaxHitDP(inta_result2.at[i, 'Demo'], int(inta_result2.at[i, 'Max_Hit_Len']))
    
# inta_result2.info(memory_usage='deep') # dataframe占用内存


# 保存CSV
# re-order columns to make sure the order is the same with RNAhybrid results
cols = ['tRF_ID', 'Transcript_ID', 'MFE', 'P_Val', 'Demo', 'Max_Hit_Len', 'Start_tRF', 'End_tRF', 'Start_Target', 'End_Target', 'Tool', 'HybridDP', 'SubseqDP', 'Max_Hit_DP']
inta_result2.to_csv(output_file, index=False, columns=cols)
print('parsed results saved to file {}!'.format(output_file))

# 整个pipeline耗时
print('Parsing IntaRNA results completed. Elapsed time: {:.2f} hours'.format((time()-begin_time)/3600.0))
