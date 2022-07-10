# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 03:07:43 2019

@author: hill103
"""

"""
多进程调用RNAhybrid程序，完成tRFs至mRNAs的绑定位置预测
多进程验证：5进程同时运行tRF 1001数据，结果与bash命令单进程结果逐行比较，结果完全一致，验证通过！
1) 删除所有数据库操作，数据先保存为CSV文件
2) 添加额外的解析函数，一次性将RNAhybrid结果转换成能存入MySQL数据库的格式
3) 一些附件信息文件改用json格式
4) 增加Tool，Gene_Evi和Site_Evi column
5) 删除所有的transcript ID和Name解析，直接输出transcript FASTA文件中的ID
输入：tRFs和transcripts序列
输出：1) trfs_infos.csv：tRFs序列信息
2) transcripts_infos.csv：transcripts序列信息
3) trf_rnahybrid_ana_infos.csv：程序运行耗时
4) rnahybrid_result.csv：解析后的binding interaction entries
"""
                                            
      
import sys, os
from shutil import rmtree
from getopt import getopt
from time import time
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
from subprocess import check_output
import pandas as pd
from datetime import datetime, timezone
from math import floor
import re


# ---------------------解析RNAhybrid图示-------------------------------
def parseDemo(demo):
    '''解析RNAhybrid的interaction图示
    返回VRNA dot-bracket notation (interaction sites only)
    返回格式为target seq&query seq，并且二者都是5'->3'方向
    '''
    
    def getBindOnly(seq, note, marker):
        '''只返回interaction区域
        marker为'('或者')'
        '''
        start = note.find(marker)
        stop = note.rfind(marker)
        return seq[start:stop+1], note[start:stop+1]
    
    # demo可分为4行，长度都相等
    # 第1和2行属于target，第3和4行属于query
    lines = demo.split('\n')
    
    target_seq = ''
    query_seq = ''
    target_note = ''
    query_note = ''
    
    # 前面9个字符为无用字符"target 5'"和"tRF 3‘"
    # 后面2个字符为无用字符"3'"和"5'"
    for i in range(9, len(lines[0])-2):
        if lines[0][i]==' ' and lines[1][i]==' ' and lines[2][i]==' ' and lines[3][i]==' ':
            # 无效列
            continue
        # 对于target或者query的1列来说，只有三种情况
        # 1.都为空白
        if lines[0][i]==' ' and lines[1][i]==' ':
            pass
        # 2.一个没有匹配上的base
        elif lines[0][i]!=' ' and lines[1][i]==' ':
            target_seq += lines[0][i]
            target_note += '.'
        # 3.一个匹配上的base
        elif lines[0][i]==' ' and lines[1][i]!=' ':
            target_seq += lines[1][i]
            target_note += '('
        else:
            raise Exception('Invalid RNAhybrid Illustration')
        
        # 同理，处理query
        # 1.都为空白
        if lines[2][i]==' ' and lines[3][i]==' ':
            pass
        # 2.一个没有匹配上的base
        elif lines[3][i]!=' ' and lines[2][i]==' ':
            query_seq += lines[3][i]
            query_note += '.'
        # 3.一个匹配上的base
        elif lines[3][i]==' ' and lines[2][i]!=' ':
            query_seq += lines[2][i]
            query_note += ')'
        else:
            raise Exception('Invalid RNAhybrid Illustration')
    
    # 只取interaction区域的sequence和notation
    target_seq, target_note = getBindOnly(target_seq, target_note, '(')
    query_seq, query_note = getBindOnly(query_seq, query_note, ')')
    
    # 用'&'连接target和query，并且query反向
    return target_seq+'&'+query_seq[::-1], target_note+'&'+query_note[::-1]


# ---------------------新Demo生成相关函数(同解析IntaRNA结果所用函数)---------------------------
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


# ---------------------其余解析用小函数-------------------------------
def getStartEnd(full_seq, sub_seq, start_index=None):
    '''计算sub_seq在full_seq中的起止坐标
    坐标从**1**开始
    '''

    if not start_index:
        start_index = 0
    # 在start_index上下游范围内检索
    # 上游范围不宜放宽，否则会找到位于其它位置的同样序列
    # 下游范围不宜过窄，否则像ts-100这样长度为47的tRF，整个匹配长度会超过100，从而报错
    start = max(0, start_index-1)
    stop = min(len(full_seq), start_index+len(sub_seq)+50)
    
    # 如果找不到sub_seq，index函数会报错
    sub_start = full_seq.index(sub_seq, start, stop)
    
    # 返回以**1**为起点的起止坐标，并转为整数形式
    return int(sub_start+1), int(sub_start+len(sub_seq)-1+1)
    
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
    for i in dataframe.index:
        if dataframe.at[i, 'Transcript_ID'] == this_enst:
            # 同一个transcript的记录，加入比较
            ind_for_check.append(i)
        else:
            # 新的transcript的开始，之前的可以进行比较了
            tmp_result = checkEntries(dataframe.loc[ind_for_check])
            if not tmp_result is None:
                if len(tmp_result) > 0:
                    '''
                    print('######Total {:d} duplicated entries for tRF "{}" and transcript "{}"######'
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
                print('######Total {:d} duplicated entries for tRF "{}" and transcript "{}"######'
                      .format(len(tmp_result), dataframe.at[tmp_result[0], 'tRF_ID'], this_enst))
                '''
                to_del += tmp_result
            
    return to_del


# ---------------------解析RNAhybrid原始结果函数-------------------------------
def getMaxHitLen(demo):
    '''给出一个匹配的示意图，找到连续匹配上的bases的最大长度
    示意图共4行，其中匹配上的序列在第2行和第3行
    '''
    tmp = demo.split('\n')
    matched = tmp[1].strip().split()
    max_len = 0
    for seq in matched:
        if max_len < len(seq):
            max_len = len(seq)
    return max_len


def getPart(pos, part_dict):
    '''根据绑定位置，以及mRNA的区域起止说明
    判断绑定点位于哪一个区域
    '''
    tmp_pos = int(pos)
    for k, v in part_dict.items():
        # 分隔起止位置
        tmp = v.split('-')
        if tmp_pos>=int(tmp[0]) and tmp_pos<=int(tmp[1]):
            return k
    
    
def parseResult(text, write_file):
    '''解析RNAhybrid本地运行返回的结果
    本地当前版本为2.1.2
    '''
    # 每一个匹配结果包含多行，分隔符为'\n\n\n'
    blocks = text.strip().split('\n\n\n')

    # 解析每一个block
    # 以分隔符'\n\n'可以将每一个block细分为3部分
    # skip CDS/UTR3/UTR5 parsing
    result = {'tRF_ID':[], 'Transcript_ID':[], 
          'MFE':[], 'P_Val':[],
          'Pos':[], 'Demo':[], 'Max_Hit_Len':[]}

    for block in blocks:
        total_parts = block.split('\n\n')
        # 第一部分包含target mRNA的ID和Length，以及miRNA的ID和Length
        # 但是可能存在一行无用信息，如target too long，其后面无空白行，且可能存在多行
        # 因此抽取信息时，使用倒数索引
        first_part = total_parts[0].split('\n')
        
        result['tRF_ID'].append(first_part[-2].split(':')[1].strip())
        # 解析mRNA编号, excluding the beginning string "target: "
        # 注意：mRNA ID中也包含':'
        # skip parsing gene Ensembl ID and other info
        result['Transcript_ID'].append(first_part[-4].split(': ')[1].strip())
        # 最后一个为空字符串，倒数2-4为'UTR5','CDS','UTR3'起止位置，但不一定3个部分都有
        # also skip parsing
       
        # 第二部分包含能量和p值
        second_part = total_parts[1].split('\n')
        # 不记录能量单位kcal/mol
        result['MFE'].append(second_part[0].split(':')[1].strip().split()[0])
        result['P_Val'].append(second_part[1].split(':')[1].strip())
        # 第三部分包含匹配位置和示意图
        tmp_index = total_parts[2].find('\n')
        result['Pos'].append(total_parts[2][:tmp_index].split()[1])
        result['Demo'].append(total_parts[2][tmp_index:].strip().replace('miRNA', 'tRF  '))
        # 计算连续匹配上的bases的最大长度
        result['Max_Hit_Len'].append(getMaxHitLen(result['Demo'][-1]))
        # 计算匹配位置属于哪一个区域
        # skipped

    # 转换成dataframe并保存
    tmp_data = pd.DataFrame.from_dict(result)
    tmp_data = tmp_data.astype({'P_Val':'float64', 'tRF_ID':str,
                 'Pos':'int64', 'MFE':'float64', 'Max_Hit_Len':'int64'})
    tmp_data.to_csv(write_file, index=False)
    return True
    
    
def rna_work(cmd):
    '''执行bash命令cmd
    '''
    start_time = time()
    print('----------------------------------')
    print('Excuted bash command: "{}"'.format(' '.join(cmd[:-1])))
    
    # check_output returns a bytes object, not a str
    # A decode process is needed first
    output = check_output(cmd[:-1]).decode('utf-8')
    '''
    # Save temporary text file
    with open(cmd[-1]+'.txt', 'w') as f:
        f.write(output)
    print('----------------------------------')
    print('Textual result file "{}" saved successfully.'.format(cmd[-1]+'.txt'))
    '''
    # Parse text string to a csv file
    parseResult(output, cmd[-1]+'.csv')
    print('----------------------------------')
    print('Parsed result file "{}" saved successfully.'.format(cmd[-1]+'.csv'))
    elapsed_time = (time()-start_time)/3600.0
    print('Elapsed time: {:.2f} hours.'.format(elapsed_time))
    return (os.path.split(cmd[-1])[-1],
            datetime.now(timezone.utc),
            elapsed_time)
    
       
def rna_analysis(target_file, query_file, output_path, n_cores, mfe=-15, mcl=6, suboptimal=1):
    '''执行tRFs与mRNA的配对
    '''
    
    # 定义最终保存文件的文件名
    # 最终保存的大CSV路径加文件名
    binding_file = os.path.join(output_path, 'rnahybrid_results.csv')
    trf_info_file = os.path.join(output_path, 'trfs_info.csv')
    run_info_file = os.path.join(output_path, 'trf_rnahybrid_ana_info.csv')
    tran_info_file = os.path.join(output_path, 'transcripts_info.csv')
    
    # 测试target mRNA文件是否为fasta格式
    with open(target_file, 'rt') as f:
        # If it's not a fasta file, any function will return false
        if not any(SeqIO.parse(f, 'fasta')):
            raise Exception('Input target file "{}" is not a valid fasta file!'.format(target_file))
    
    # 测试query tRFs文件是否为fasta格式
    with open(query_file, 'rt') as f:
        # If it's not a fasta file, any function will return false
        if not any(SeqIO.parse(f, 'fasta')):
            raise Exception('Input query file "{}" is not a valid fasta file!'.format(query_file))
        
    # 在输出路径下建立temporary directory
    directory = os.path.join(output_path, 'RNA_tmp_files')
    if os.path.isdir(directory):
        # 删除文件夹后，再新建
        rmtree(directory)
        print('WARNING: directory "{}" removed!'.format(directory))
        
    os.mkdir(directory)
    print('Temporary directory "{}" created!'.format(directory))
    
    # 将每一条query sequence保存成一个单独的fasta文件，以实现并行计算
    output_list = []
    id_dict = {}
    count = 0
    cmds = [] # 需要执行的bash命令
    tRF_info = []
    with open(query_file, 'rt') as f:
        for record in SeqIO.parse(f, 'fasta'):
            tmp_query_file = os.path.join(directory, 'tRF_{:d}.fasta'.format(count))
            tmp_output_file = os.path.join(directory, 'tRF_{:d}_result'.format(count))
            id_dict['tRF_{:d}_result'.format(count)] = record.id
            output_list.append(tmp_output_file+'.csv')
            SeqIO.write(record, tmp_query_file, 'fasta')
            count += 1
            cmds.append(['/app/RNAhybrid', '-q', tmp_query_file, '-t', target_file,
                             '-b', str(suboptimal), '-e', str(mfe), '-m', str(150000),
                             '-n', str(70), '-s', '3utr_human', tmp_output_file])
            # tRF_ID是unique的
            tRF_info.append({'tRF_ID':record.id,
                             'tRF_Seq':str(record.seq).strip(),
                             'tRF_Length': len(record.seq.strip())})
    print('Total {:d} tRF sequences saved in temporary directory "{}".'.format(count, directory))
    
    # tRF信息存入CSV文件
    pd.DataFrame(tRF_info).to_csv(trf_info_file, index=False)
    # 特别建立一个tRFs序列的dict，方便查询，并且'T'变成'U'
    trf_seq = {}
    for item in tRF_info:
        trf_seq[item['tRF_ID']] = item['tRF_Seq'].replace('T', 'U')
    # check whether all tRF IDs are unique
    assert(len(trf_seq) == count)
        
    del tRF_info
    
    # Python调用bash命令，执行并行计算
    print('Total number of CPUs: {:d}.'.format(cpu_count()))
    print('{:d} CPUs will be used for RNAhybrid.'.format(n_cores))
    pool = Pool(n_cores)
    results = pool.map(rna_work, cmds)
    # 关闭线程池，等待工作结束
    pool.close()
    pool.join()
    
    
    # 程序运行耗时表存入CSV文件
    '''
    pd.DataFrame([{'tRF_ID':id_dict[a], 'UTC_Time':b, 'Used_Time_Hours':c}
        for a,b,c in results]).to_csv(run_info_file, index=False)
    '''

    # 准备附加信息
    # 1.tRF信息dict tRF_infos，用tRF ID索引，已存入文件
    
    # 2.gene symbol信息dict gene_symbol，用gene ensembl id索引
    # skipped
    
    # 3.transcript type和name信息rna_infos，用transcript ensembl id索引
    # skipped
    
    # 4.transcript 序列信息rna_seq，用transcript ensembl id索引
    # skip transcript parsing
    rna_seq = []
    with open(target_file, 'rt') as f:
        for record in SeqIO.parse(f, 'fasta'):
            # skip transcript ID parsing
            rna_seq.append({'Trans_ID': record.id,
                            'Trans_Seq': str(record.seq).strip(),
                            'Trans_Length': len(record.seq.strip())})

    # transcript信息写入CSV文件
    pd.DataFrame(rna_seq).to_csv(tran_info_file, index=False)
    # 特别建立一个transcript序列的dict，方便查询，并且'T'变成'U'
    tran_seq = {}
    for item in rna_seq:
        tran_seq[item['Trans_ID']] = item['Trans_Seq'].replace('T', 'U')
    # check whether all tRF IDs are unique
    assert(len(tran_seq) == len(rna_seq))
    
    del rna_seq
    
    # 将所有csv文件进一步解析，并合并成一个大CSV文件
    print('Combining results of each tRF...')
    for ind, one_file in enumerate(output_list):
        start_time = time()
        print('Processing file "{}"...'.format(one_file))
        data = pd.read_csv(one_file)
            
        # 抽取匹配长度大于等于6的序列
        print('Total {:,} entries'.format(data.shape[0]))
        data = data[data['Max_Hit_Len']>=mcl]
        print('After exclude entries with max_hit_len<{:d}, remaining {:,} entries'.format(mcl, data.shape[0]))
        
        # 进一步解析结果
        for i in data.index:
           # 从RNAhybrid的Demo中抽取interaction图示
           data.at[i, 'subseqDP'], data.at[i, 'hybridDP'] = parseDemo(data.at[i, 'Demo'])
           # 确定interaction在tRF和transcript上的起止坐标
           # 坐标以**1**为起点
           # RNAhybrid的结果，pos从0开始
           # 注意tRF的interaction序列无需反转（本来就是5'->3'方向）
           data.at[i, 'Start_Target'], data.at[i, 'End_Target'] = getStartEnd(
                tran_seq[data.at[i, 'Transcript_ID']], data.at[i, 'subseqDP'].split('&')[0],
                data.at[i, 'Pos'])
           data.at[i, 'Start_tRF'], data.at[i, 'End_tRF'] = getStartEnd(
                trf_seq[data.at[i, 'tRF_ID']], data.at[i, 'subseqDP'].split('&')[1])
           # 生成新demo
           data.at[i, 'Demo'] = getDemo(tran_seq[data.at[i, 'Transcript_ID']],
               int(data.at[i, 'Start_Target']), int(data.at[i, 'End_Target']),
               trf_seq[data.at[i, 'tRF_ID']],
               int(data.at[i, 'Start_tRF']), int(data.at[i, 'End_tRF']),
               data.at[i, 'subseqDP'], data.at[i, 'hybridDP'])
           data.at[i, 'Max_Hit_DP'] = getMaxHitDP(data.at[i, 'Demo'], int(data.at[i, 'Max_Hit_Len']))
    
        # 确认重复entries
        # 根据时间测试结果，该功能通常条件下需要花费1小时，需要进行优化，避免重复进行dataframe的column操作
        # 优化后耗时通常条件下为2分钟
        
        # 对所有entries进行排序
        data.sort_values(['Transcript_ID', 'MFE', 'Max_Hit_Len'],
                     ascending=[True, True, False], inplace=True)
    
        to_del = checkDuplicate(data)
        
        print('detected {:d} duplicated interactions'.format(len(to_del)))
        
        # 删除重复entries by index
        if len(to_del) > 0:
            data.drop(to_del, inplace=True)
    
        # 准备保存成大CSV文件
        data.drop(columns=['Pos'], inplace=True)
        data.rename(columns={'hybridDP': 'HybridDP', 'subseqDP': 'SubseqDP'}, inplace=True)
        data['Tool'] = 'RNAhybrid'
        
        # recorde all columns for saving
        # update: do not save p value 
        cols = ['tRF_ID', 'Transcript_ID', 'MFE', 'Demo', 'Max_Hit_Len', 'Start_tRF', 'End_tRF', 'Start_Target', 'End_Target', 'Tool', 'HybridDP', 'SubseqDP', 'Max_Hit_DP']
        
        # if file does not exist write header
        if ind==0:
            data.to_csv(binding_file, header=True, columns=cols, index=False)
        else: # else it exists so append without writing the header
            data.to_csv(binding_file, mode='a', header=False, columns=cols, index=False)
        
        print('Elapsed time: {:.2f} hours.'.format((time()-start_time)/3600.0))

    print('All csv results inserted into big file "{}"'.format(binding_file))
    
    
    # 删除临时文件夹
    rmtree(directory)
    print('WARNING: temporary directory "{}" removed!'.format(directory))
    
    
    return True


#############主函数#####################################################################
def usage():
    '''对主函数进行简介
    '''
    print('''
python callRNAhybrid.py [option][value]...
    -h or --help        print this help messages
    -t or --target      mRNA fasta file for input，with absolute or relative path
    -q or --query       tRFs fasta file for input, with absolute or relative path
    -o or --outputpath  absolute or relative path for result files. If ignored, the current path will be used
    -n or --n_cores     number of CPU cores used for parallel computation. Default value is 1
    -e or --MFE         free energy threshold, used for RNAhybrid `-e` option. Default value is -15
    -m or --MCL         threshold of maximum complementary length, and interactions with maximum complementary length less than it are filtered out. Default value is 6
    -b or --suboptimal  reported number of interaction sites on each transcript, used for RNAhybrid `-b` option. Default value is 1
''')

# 如果没有任何参数，显示提示信息，并退出
if len(sys.argv) == 1:
    print('-h or --help for detail')
    sys.exit(1)


# 定义命令行参数
# 短选项名后的冒号(:)表示该选项必须有附加的参数
# 长选项名后的等号(=)表示该选项必须有附加的参数。
shortargs = 'ht:q:o:n:e:m:b:'
longargs = ['help', 'target=', 'query=', 'outputpath=', 'n_cores=', 'MFE=', 'MCL=', 'suboptimal=']

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
paramdict = {'target_file':None, 'query_file':None, 'output_path':os.getcwd(),
             'n_cores':1, 'MFE':-15, 'MCL':6, 'suboptimal':1}
 
for opt,val in opts:
    if opt in ('-h', '--help'):
        usage()
        sys.exit(1)
        
    if opt in ('-t', '--target'):
        if not os.path.isfile(val):
            # 输入不是一个确实存在的文件名
            raise Exception('Invalid input target file!')
        # 采用realpath函数，获得真实绝对路径
        paramdict['target_file'] = os.path.realpath(val)
        continue
    
    if opt in ('-q', '--query'):
        if not os.path.isfile(val):
            # 输入不是一个确实存在的文件名
            raise Exception('Invalid input query file!')
        # 采用realpath函数，获得真实绝对路径
        paramdict['query_file'] = os.path.realpath(val)
        continue
    
    if opt in ('-o', '--outputpath'):
        if not os.path.isdir(val):
            # 输入不是一个确实存在的文件夹
            raise Exception('Invalid or not existed directory!')
        # 采用realpath函数，获得真实绝对路径
        paramdict['output_path'] = os.path.realpath(val)
        continue
    
    if opt in ('-n', '--n_cores'):
        paramdict['n_cores'] = int(val)
        continue
    
    if opt in ('-e', '--MFE'):
        paramdict['MFE'] = float(val)
        continue
    
    if opt in ('-m', '--MCL'):
        paramdict['MCL'] = int(val)
        continue
    
    if opt in ('-b', '--suboptimal'):
        paramdict['suboptimal'] = int(val)
        continue

# 检查参数是否齐全
for k,v in paramdict.items():
    if v is None:
        raise Exception('Option "{}" is missing!'.format(k))

print('RNAhybrid running options:')
for k,v in paramdict.items():
    print('{}: {}'.format(k, v))

# 调用分析函数
start_time = time()
rna_analysis(paramdict['target_file'], paramdict['query_file'], paramdict['output_path'],
             paramdict['n_cores'], paramdict['MFE'], paramdict['MCL'], paramdict['suboptimal'])
print('RNAhybrid analysis completed. Elapsed time: {:.2f} hours.'.format((time()-start_time)/3600.0))
