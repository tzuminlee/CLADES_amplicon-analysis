from __future__ import division
from collections import Counter
import numpy as np
from pandas import Series, DataFrame
import pandas as pd
import copy
import sys
import cPickle as pickle
from tqdm import tqdm

def save_object(obj, filename, n = 6000):
	sys.setrecursionlimit(n)
	with open(filename, 'wb') as output:
		pickle.dump(obj, output, -1)

def retrieve_object(filename):
	with open(filename, 'rb') as input:
		obj = pickle.load(input)
	return obj

def LCS_score(X, Y):
	n, m = len(X), len(Y)
	L = [[0] * (m+1) for k in range(n+1)]
	for j in range(n):
		for k in range(m):
			if X[j] == Y[k]:
				L[j+1][k+1] = L[j][k] + 1
			else:
				L[j+1][k+1] = max(L[j][k+1], L[j+1][k])
	arr = np.array(L)
	return arr.max(axis=1).max()

def codon_reverse(seq):
	seq_list = list(seq)
	for i, j in enumerate(seq_list):
		if j == 'A': seq_list[i] = 'T'
		if j == 'T': seq_list[i] = 'A'
		if j == 'G': seq_list[i] = 'C'
		if j == 'C': seq_list[i] = 'G'
	seq_list.reverse()
	seq_reverse = ''.join(seq_list)
	return seq_reverse

def string_reverse(string):
	string_list = list(string)
	string_list.reverse()
	return ''.join(string_list)

def read_score_mean(score_string):
	total_score = 0
	for i in score_string:
		total_score += ord(i)
	return total_score/len(score_string)

def FASTQ_qc1(R1_filepath, R2_filepath):
	R1_seq, R1_score, R2_seq, R2_score = [], [], [], []
	with open(R1_filepath, 'rb') as test:
		first = 1
		last = 1000000
		count = 0
		for line in tqdm(test):
			count += 1
			if count > last*4: break
			if count >  (first-1)*4:
				if count % 4 == 2:
					R1_seq.append(line)
				if count % 4 == 0:
					R1_score.append(line)
	with open(R2_filepath, 'rb') as test:
		first = 1
		last = 1000000
		count = 0
		for line in tqdm(test):
			count += 1
			if count > last*4: break
			if count >  (first-1)*4:
				if count % 4 == 2:
					R2_seq.append(line)
				if count % 4 == 0:
					R2_score.append(line)
	R2_seq_rev = []
	for i in R2_seq:
		R2_seq_rev.append(codon_reverse(i))
	R2_score_rev = [string_reverse(i) for i in R2_score]
	R1_score_mean = [read_score_mean(i) for i in R1_score]
	arr1 = np.array(R1_score_mean)
	R1_mean, R1_std = arr1.mean(), arr1.std()
	R1_score_th = arr1.mean() - arr1.std()*4
	print
	print R1_mean, R1_std, R1_score_th
	R1_fail = 0
	R1_pass = 0
	for i in R1_score_mean:
		if i < R1_score_th:
			R1_fail += 1
		else: R1_pass += 1
	print len(R1_seq), R1_fail/len(R1_seq), R1_pass/len(R1_seq), (R1_fail+R1_pass)/len(R1_seq)
	R2_score_mean = [read_score_mean(i) for i in R2_score_rev]
	arr2 = np.array(R2_score_mean)
	R2_mean, R2_std = arr2.mean(), arr2.std()
	R2_score_th = arr2.mean() - arr2.std()*4
	print R2_mean, R2_std, R2_score_th
	R2_fail = 0
	R2_pass = 0
	for i in R2_score_mean:
		if i < R2_score_th:
			R2_fail += 1
		else: R2_pass += 1
	print len(R2_seq_rev), R2_fail/len(R2_seq_rev), R2_pass/len(R2_seq_rev), (R2_fail+R2_pass)/len(R2_seq_rev)
	print
	R1_select1, R1_select1_score, R2_select1, R2_select1_score = [], [], [], []
	for i, j in enumerate(R1_seq):
		k = R2_seq_rev[i]
		if 'N' not in j and 'N' not in k and len(j) > 30 and len(k) > 30:
			if R1_score_mean[i]>R1_score_th and R2_score_mean[i]>R2_score_th:
				R1_select1.append(j)
				R1_select1_score.append(R1_score[i])
				R2_select1.append(k)
				R2_select1_score.append(R2_score_rev[i])
	print len(R1_seq), len(R1_select1), len(R1_select1)/len(R1_seq)
	print
	return R1_select1, R1_select1_score, R2_select1, R2_select1_score	

def stitch_reads5(R1, R1_score, R2, R2_score):
	assembled = []
	for i in tqdm(range(len(R1))):
		r1, r1_score, r2, r2_score = R1[i], R1_score[i], R2[i], R2_score[i]
		query = r2[:len(r2)//5]
		arr0 = DR_matrix(r1, query)[0]
		arr = np.array(arr0)
		if arr.max() > (len(r2)//5)*0.9:
			idx = np.unravel_index(arr.argmax(), arr.shape)
			answer = r1[:idx[0]] + r2[idx[1]:]
			assembled.append(answer)
		else:
			query = r1[-(len(r1)//5):]
			arr0 = DR_matrix(r2, query)[0]
			arr = np.array(arr0)
			if arr.max() > (len(r1)//5)*0.9:
				idx = np.unravel_index(arr.argmax(), arr.shape)
				answer = r1[:len(r1)-((len(r1)//5)-idx[1])+1] + r2[idx[0]:]
				assembled.append(answer)
	print len(R1), len(assembled), len(assembled)/len(R1)
	print
	return assembled

def DR_matrix(S1, S2):
	n, m = len(S1), len(S2)
	arr = [[0]*m for i in range(n)]
	for i in range(n):
		for j in range(m):
			if S1[i] == S2[j]:
				arr[i][j] = 1
	for i in range(n):
		for j in range(m):		
			if arr[i][j] == 1:
				if i > 0 and j > 0:
					arr[i][j] = arr[i-1][j-1] + 1	
	arr_2 = copy.deepcopy(arr)
	for i in range(n):
		j = 0
		while j < m:
			a = 0
			if arr_2[i][j] == 1:
				while i+1 < n and j+1 < m:
					if arr_2[i+1][j+1] - arr_2[i][j] == 1:
						a += 1
						last_value = arr_2[i+1][j+1]
						i += 1
						j += 1
					else: break
			if a > 0:
				for k in range(a):
					i -= 1
					j -= 1
					arr_2[i][j] = last_value
			j += 1
	return arr, arr_2

def align1(Ori_seq, LCS_seq, index_arr, seg_list=[], backward=False):
	a, score_table = DR_matrix(Ori_seq, LCS_seq)
	arr = np.array(score_table)
	largest = 0
	if backward:
		for i in range(arr.shape[0]-1, -1, -1):
			for j in range(arr.shape[1]-1, -1, -1):
				if arr[i][j] > largest:
					largest = arr[i][j]
					index_start = index_arr[i-largest+1][j-largest+1]
					index_end = index_arr[i][j]
					largest_start = (i-largest+1, j-largest+1)
					largest_end = (i, j)
	else:
		for i in range(arr.shape[0]):
			for j in range(arr.shape[1]):
				if arr[i][j] > largest:
					largest = arr[i][j]
					index_start = index_arr[i][j]
					index_end = index_arr[i+largest-1][j+largest-1]
					largest_start = (i, j)
					largest_end = (i+largest-1, j+largest-1)
	if largest < 1:
		seg_list.append([])
	elif largest == arr.shape[1]:
		seg_list.append([index_start, index_end])
	else:
		if largest_start[1] > 0:
			align1(Ori_seq[:largest_start[0]], LCS_seq[:largest_start[1]], index_arr[:largest_start[0], :largest_start[1]], seg_list, backward=True)
		seg_list.append([index_start, index_end])
		if largest_end[1] < arr.shape[1] - 1:
			align1(Ori_seq[largest_end[0]+1:], LCS_seq[largest_end[1]+1:], index_arr[largest_end[0]+1:, largest_end[1]+1:], seg_list, backward=False)
	return seg_list

def find_contiguous3(seq, query):
	n, m = len(seq), len(query)
	index_arr = np.arange(0, n*m, 1).reshape(n, m)
	seg_list = align1(seq, query, index_arr, seg_list=[])
	seg_index = []
	empty = False
	for i in seg_list:
		if i:
			seg_index.append([(int(i[0]//m), i[0]%m), (int(i[1]//m), i[1]%m)])
		else:
			empty = True
	seg_index_copy = [i for i in seg_index]
	seg_index_rev = []
	x = 0
	if len(seg_index) < 2:
		seg_index_rev = seg_index
	else:
		while x+1 < len(seg_index):
			a, b = seg_index[x][1][0], seg_index[x+1][0][0]
			j, k = seg_index[x][1][1], seg_index[x+1][0][1]
			c = 0
			while b-a+c > 1:
				if seq[a+1+c] == seq[b+c] and query[j+1+c] == query[k+c]:
					c += 1
				else: break
			if c > 0:
				seg_index_rev.append([seg_index[x][0], (seg_index[x][1][0]+c, seg_index[x][1][1]+c)])
				seg_index[x+1] = [(seg_index[x+1][0][0]+c, seg_index[x+1][0][1]+c), seg_index[x+1][1]]
				if x + 2 == len(seg_index):
					seg_index_rev.append(seg_index[-1])
			else:
				seg_index_rev.append(seg_index[x])
				if x + 2 == len(seg_index):
					seg_index_rev.append(seg_index[-1])		
			x += 1	
	return seg_index_copy, seg_index_rev, empty	

def annotating_seq(seq, index, template):
	seq_list = []
	i = 0
	if index[0][0][0] > 0:
		seq_list.append(str(index[0][0][0]))
	if index[0][0][1] > 0:
		seq_list.append(seq[:index[0][0][1]])
	if i+1 == len(index):
		seq_list.append(seq[index[0][0][1]:index[0][1][1]+1])
		if index[0][1][1]+1 < len(seq):
			seq_list.append(seq[index[0][1][1]+1:])
	else:
		while i+1 < len(index):
			j, k = index[i], index[i+1]
			seq_list.append(seq[j[0][1]:j[1][1]+1])
			gap = k[0][0] - j[1][0] - 1
			if gap > 0:
				seq_list.append(str(gap))
			extra = k[0][1] - j[1][1] - 1
			if extra > 0:
				seq_list.append(seq[j[1][1]+1:k[0][1]])
			if i+2 == len(index):
				seq_list.append(seq[k[0][1]:k[1][1]+1])
				if k[1][1]+1 < len(seq):
					seq_list.append(seq[k[1][1]+1:])
				if k[1][0]+1 < len(template):
					seq_list.append(str(len(template)-(k[1][0]+1)))
			i += 1
	seq_annotated = '-'.join(seq_list)
	return seq_annotated

def remove_minor(homology_index, minor_th):
	answer = []
	if homology_index:
		for i in homology_index:
			if i[1][0] - i[0][0] > minor_th:
				answer.append(i)
	return answer

def indel_index5(index_list, read, seq):
	deletion, insertion, indel, mismatch = [], [], [], []
	n = 0
	while n+1 < len(index_list):
		a, b, c, d = index_list[n][1][0], index_list[n][1][1], index_list[n+1][0][0], index_list[n+1][0][1]
		if c - a > 1 and d == b+1:
			deletion.append(((a+1, c-1),))
		if d - b > 1 and c == a+1:
			m = 0
			while seq[c+m] == read[b+1+m]:
				m += 1
			insertion.append(((a+m, c+m), d-b-1, (read[b+1+m:d+m])))
		if c - a > 1 and d - b > 1:
			if c-a != d-b:
				m = 0
				while seq[c+m] == read[b+1+m]:
					m += 1
				indel.append(((a+1+m, c-1+m), (b+1+m, d-1+m)))
			else:
				mismatch.append(((a+1, c-1), (b+1, d-1)))
		n += 1
	return deletion, insertion, indel, mismatch

def find_SSA(del_index, seq):
	a, b = del_index[0], del_index[1]
	n, m = 0, 0
	while a+n <= b and seq[a+n] == seq[b+n+1]:
		n += 1
	while b-m >= a and seq[b-m] == seq[a-m-1]:
		m += 1
	return tuple((n, m))

def find_UMI(read, seq1, seq2):
	query = None
	if seq1 in read[:80] and seq2 in read[:80]:
		query = read[:80]
	if seq1 in read[-80:] and seq2 in read[-80:]:
		query = read[-80:]
	if query:
		index1, index2, empty1 = find_contiguous3(seq1, query)
		index3, index4, empty2 = find_contiguous3(seq2, query)	
		return query[index1[-1][-1][-1]+1:index3[0][0][1]]
	else: return None

def find_UMI2(read, seq):	
	query = None
	if seq in read[:50]:
		query = read[:50]
	if seq in read[-50:]:
		query = read[-50:]
	if query:
		query_list = list(query)
		query = ''.join(a for a in query_list if a.isalpha())		
		index1, index2, empty1 = find_contiguous3(seq, query)
		return query[:index1[0][0][1]]

def sort_UMI(UMI_dict):
	UMI_list = []
	for k, v in UMI_dict.items():
		UMI_list.append((k, v))
	merge_sort_UMIs(UMI_list)
	UMI_sorted = [i[0] for i in UMI_list]
	return UMI_sorted

def merge_sort_UMIs(S):
	n = len(S)
	if n < 2:
		return
	mid = n // 2
	S1 = S[0:mid]
	S2 = S[mid:n]
	merge_sort_UMIs(S1)
	merge_sort_UMIs(S2)
	merge_UMIs(S1, S2, S)

def merge_UMIs(S1, S2, S):
	i = j = 0
	while i + j < len(S):
		if j == len(S2) or (i < len(S1) and len(S1[i][1]) > len(S2[j][1])):
			S[i+j] = S1[i]
			i += 1
		else:
			S[i+j] = S2[j]
			j += 1


def main1(sample_ID, date, folder_path):
	R1_filepath = folder_path + sample_ID + '_R1_001.fastq'
	R2_filepath = folder_path + sample_ID + '_R2_001.fastq'	
	R1_seq, R1_score, R2_seq, R2_score = FASTQ_qc1(R1_filepath, R2_filepath)
	assembled = stitch_reads5(R1_seq, R1_score, R2_seq, R2_score)
	save_object(assembled, folder_path + sample_ID + '_assembled_' + date + '.pkl')
	counted = Counter(assembled)
	print counted.most_common(3)
	print

def main2(sample_ID, date, folder_path):
	assembled = retrieve_object(folder_path + sample_ID + '_assembled_' + date + '.pkl')
	read_dict1 = {}
	for i in tqdm(assembled):
		i_list = list(i)
		i = ''.join(a for a in i_list if a.isalpha())
		i_UMI = find_UMI2(i, 'CGCCAAGCAGAGAGGGCG')
		if i_UMI:
			if i_UMI not in read_dict1.keys():
				read_dict1[i_UMI] = [i]
			else: read_dict1[i_UMI].append(i)
	UMI_sorted = sort_UMI(read_dict1)
	read_dict2 = {}
	for k, v in read_dict1.items():
		counted = Counter(v)
		ranked = counted.most_common()
		read_dict2[k] = ranked
	read_list = []
	length, count = [], []
	spacer2_seq = 'CATCCATACAGTACCA'
	spacer2 = []
	UMI_list, UMI_length, UMI_count = [], [], []
	for i in UMI_sorted:
		m = 0
		for s, c in read_dict2[i]:
			if c > 0:
				m += 1
		n = 1
		for s, c in read_dict2[i]:
			if c > 0:			
				if spacer2_seq in s:
					spacer2.append('yes')
				else: spacer2.append('no')
				read_list.append(s)
				length.append(len(s))
				count.append(c)
				UMI_list.append((i, n))
				UMI_length.append(len(i))
				UMI_count.append(m)
				n += 1
	temp1, temp2, temp3, temp4, temp5, temp6, temp7 = Series(read_list), Series(length), Series(count), Series(UMI_list), Series(UMI_length), Series(UMI_count), Series(spacer2)
	df0 = pd.concat([temp1, temp2, temp3, temp4, temp5, temp6, temp7], axis=1)
	df0.columns = ['sequence', 'length', 'count', 'UMI', 'UMI_length', 'UMI_count', 'spacer2?']
	writer = pd.ExcelWriter(folder_path + sample_ID + '_assembled_' + date + '.xlsx')
	df0.to_excel(writer, 'assembled')
	df1 = df0[df0['UMI_count']>1]
	df1.to_excel(writer, 'multi_seq')
	df2 = df0[df0['UMI_count']==1]
	df2.to_excel(writer, 'uni_seq')
	save_object(read_dict2, folder_path + sample_ID + '_assembled_dict2_' + date + '.pkl')
	save_object(df0, folder_path + sample_ID + '_assembled_df0_' + date + '.pkl')
	spacer2_list = []
	length2, count2 = [], []
	spacer2_seq = 'CATCCATACAGTACCA'
	UMI_list2, UMI_length2 = [], []
	for i in UMI_sorted:
		for s, c in read_dict2[i]:
			if c > 0 and spacer2_seq in s:
				spacer2_list.append(s)
				length2.append(len(s))
				count2.append(c)
				UMI_list2.append(i)
				UMI_length2.append(len(i))
				break
	temp1, temp2, temp3, temp4, temp5 = Series(spacer2_list), Series(length2), Series(count2), Series(UMI_list2), Series(UMI_length2)
	df_spacer2 = pd.concat([temp1, temp2, temp3, temp4, temp5], axis=1)
	df_spacer2.columns = ['sequence', 'length', 'count', 'UMI', 'UMI_length']
	df_spacer2.to_excel(writer, 'spacer2')
	save_object(df_spacer2, folder_path + sample_ID + '_df_spacer2_' + date + '.pkl')

def main3(sample_ID, date, folder_path, template_seq):
	df1 = retrieve_object(folder_path + sample_ID + '_assembled_df0_' + date + '.pkl')
	df2 = df1[df1['spacer2?'] == 'yes']
	writer = pd.ExcelWriter(folder_path + sample_ID + '_assembled_df2_' + date + '.xlsx')
	queries = df2['sequence'].values
	contiguous1, contiguous2 = [], []
	with_gap = []
	for query in tqdm(queries):
		contiguous_index1, contiguous_index2, empty = find_contiguous3(template_seq, query)
		contiguous1.append(contiguous_index1)
		contiguous2.append(contiguous_index2)
		with_gap.append(empty)
	contiguous3 = []
	for i in contiguous2:
		contiguous3.append(len(i))
	contiguous4 = []
	for i in contiguous2:
		contiguous4.append(remove_minor(i, minor_th=4))
	df2.loc[:, 'contiguous1'] = contiguous1
	df2.loc[:, 'contiguous2'] = contiguous2
	df2.loc[:, 'with_gap'] = with_gap
	df2.loc[:, 'seg#'] = contiguous3
	df2.loc[:, '>4'] = contiguous4
	df2.to_excel(writer, 'uni_seq')	
	save_object(df2, folder_path + sample_ID + '_assembled_df2_' + date + '.pkl')
	print 'end'	

def main4(sample_ID, date, folder_path, template_seq):
	df2 = retrieve_object(folder_path + sample_ID + '_assembled_df2_' + date + '.pkl')
	contiguous = df2['>4'].values
	print contiguous[:5]
	reads, counts = df2['sequence'].values, df2['count'].values
	deletions, insertions, indels, mismatches = [], [], [], []
	for i, j in zip(contiguous, reads):
		a, b, c, d = indel_index5(i, j, template_seq)
		deletions.append(a)
		insertions.append(b)
		indels.append(c)
		mismatches.append(d)
	df2.insert(8, column='deletion', value=deletions)
	df2.insert(9, column='insertion', value=insertions)
	df2.insert(10, column='indel', value=indels)
	df2.insert(11, column='mismatch', value=mismatches)
	indel_count = []
	for i in df2.index:
		a, b, c = len(df2.ix[i, 'deletion']), len(df2.ix[i, 'insertion']), len(df2.ix[i, 'indel'])
		indel_count.append(a+b+c)
	df2.insert(12, column='indel_count', value=indel_count)
	seq_annotated = []
	for i, j in zip(df2['sequence'].values, df2['>4'].values):
		seq_annotated.append(annotating_seq(i, j, template_seq))
	df2.insert(13, column='annotated_seq', value=seq_annotated)
	writer = pd.ExcelWriter(folder_path + sample_ID + '_assembled_df0_indels_' + date + '.xlsx')
	df2.to_excel(writer, 'indels')
	save_object(df2, folder_path + sample_ID + '_assembled_df0_indels_' + date + '.pkl')

def main5(sample_ID, date, folder_path, template_seq):
	df2 = retrieve_object(folder_path + sample_ID + '_assembled_df0_indels_' + date + '.pkl')
	df2 = df2.sort_values(by='count', ascending=False)
	counts, length, seq_annotated, indel_count = df2['count'].values, df2['length'].values, df2['annotated_seq'].values, df2['indel_count'].values
	deletions, insertions, indels, mismatches = df2['deletion'].values, df2['insertion'].values, df2['indel'].values, df2['mismatch'].values
	writer = pd.ExcelWriter(folder_path + sample_ID + '_assembled_indels_' + date + '.xlsx')
	indel_dfs = {}
	wt_reads, multi = {}, {}
	total = 0
	for i, c in enumerate(counts):
		total += c
		if indel_count[i] == 0:
			if str(length[i]) not in wt_reads.keys():
				wt_reads[str(length[i])] = [[seq_annotated[i]], [c]]
			else:
				wt_reads[str(length[i])][0].append(seq_annotated[i])
				wt_reads[str(length[i])][1].append(c)
		if indel_count[i] > 1:
			a = (tuple(deletions[i]), tuple(insertions[i]), tuple(indels[i]))
			if a not in multi.keys():
				multi[a] = [[seq_annotated[i]], [c]]
			else:
				multi[a][0].append(seq_annotated[i])
				multi[a][1].append(c)
	if wt_reads:
		wt_reads2 = {}
		for k, v in wt_reads.items():
			count_sum = sum(v[1])
			if len(v[0]) > 2:
				first, second, third = v[0][0], v[0][1], v[0][2]
				first_count, second_count, third_count = v[1][0], v[1][1], v[1][2]
				wt_reads2[k] = [count_sum, first, first_count, second, second_count, third, third_count]
			elif len(v[0]) == 2:
				first, second, third = v[0][0], v[0][1], str()
				first_count, second_count, third_count = v[1][0], v[1][1], 0
				wt_reads2[k] = [count_sum, first, first_count, second, second_count, third, third_count]
			else:
				first, second, third = v[0][0], str(), str()
				first_count, second_count, third_count = v[1][0], 0, 0
				wt_reads2[k] = [count_sum, first, first_count, second, second_count, third, third_count]					
		df2_wt = DataFrame.from_dict(wt_reads2, orient='index')
		df2_wt.columns = ['total_count', 'first', 'first_count', 'second', 'second_count', 'third', 'third_count']
		percent = [v/total for v in df2_wt['total_count'].values]
		df2_wt.insert(1, column='frequency', value=percent)
		df2_wt_sorted = df2_wt.sort_values(by='total_count', ascending=False)	
		df2_wt_sorted.to_excel(writer, 'wt_read')
		indel_dfs['wt_read'] = df2_wt_sorted
	if multi:
		multi2 = {}
		for k, v in multi.items():
			count_sum = sum(v[1])
			if len(v[0]) > 2:
				first, second, third = v[0][0], v[0][1], v[0][2]
				first_count, second_count, third_count = v[1][0], v[1][1], v[1][2]
				multi2[k] = [count_sum, first, first_count, second, second_count, third, third_count]
			elif len(v[0]) == 2:
				first, second, third = v[0][0], v[0][1], str()
				first_count, second_count, third_count = v[1][0], v[1][1], 0
				multi2[k] = [count_sum, first, first_count, second, second_count, third, third_count]
			else:
				first, second, third = v[0][0], str(), str()
				first_count, second_count, third_count = v[1][0], 0, 0
				multi2[k] = [count_sum, first, first_count, second, second_count, third, third_count]					
		df2_multi = DataFrame.from_dict(multi2, orient='index')
		df2_multi.columns = ['total_count', 'first', 'first_count', 'second', 'second_count', 'third', 'third_count']
		percent = [v/total for v in df2_multi['total_count'].values]
		df2_multi.insert(1, column='frequency', value=percent)
		df2_multi_sorted = df2_multi.sort_values(by='total_count', ascending=False)
		indel_dfs['multi_indel'] = df2_multi_sorted

	for items, label in zip([deletions, insertions, indels, mismatches], ['deletion', 'insertion', 'indel', 'mismatch']):
		items_dict = {}
		a = 0
		total = 0
		for i, n in tqdm(zip(items, counts)):
			total += n
			if i:
				for k in i:
					if k not in items_dict.keys():
						items_dict[k] = [[seq_annotated[a]], [n]]
					else:
						items_dict[k][0].append(seq_annotated[a])
						items_dict[k][1].append(n)
			a += 1
		items_dict2 = {}
		for k, v in items_dict.items():
			count_sum = sum(v[1])
			if len(v[0]) > 2:
				first, second, third = v[0][0], v[0][1], v[0][2]
				first_count, second_count, third_count = v[1][0], v[1][1], v[1][2]
				items_dict2[k] = [count_sum, first, first_count, second, second_count, third, third_count]
			elif len(v[0]) == 2:
				first, second, third = v[0][0], v[0][1], str()
				first_count, second_count, third_count = v[1][0], v[1][1], 0
				items_dict2[k] = [count_sum, first, first_count, second, second_count, third, third_count]
			else:
				first, second, third = v[0][0], str(), str()
				first_count, second_count, third_count = v[1][0], 0, 0
				items_dict2[k] = [count_sum, first, first_count, second, second_count, third, third_count]
		if items_dict2:
			df2_items = DataFrame.from_dict(items_dict2, orient='index')
			df2_items.columns = ['total_count', 'first', 'first_count', 'second', 'second_count', 'third', 'third_count']
			percent = [v/total for v in df2_items['total_count'].values]
			df2_items.insert(1, column='frequency', value=percent)
			df2_items_sorted = df2_items.sort_values(by='total_count', ascending=False)	
			if label == 'deletion':
				size, SSA = [], []
				for i in df2_items_sorted.index:
					size.append(i[0][1]-i[0][0]+1)
					SSA.append(find_SSA(i[0], template_seq))
				df2_items_sorted.insert(0, column='SSA', value=SSA)
				df2_items_sorted.insert(1, column='size', value=size)
			if label == 'indel':
				size, SSA = [], []
				for i in df2_items_sorted.index:
					size.append(i[0][1]-i[0][0]+1)
					SSA.append(find_SSA(i[0], template_seq))
				df2_items_sorted.insert(0, column='SSA', value=SSA)
				df2_items_sorted.insert(1, column='del_size', value=size)	
			df2_items_sorted.to_excel(writer, label)
			indel_dfs[label] = df2_items_sorted
	if multi:
		df2_multi_sorted.to_excel(writer, 'multi_indel')
	save_object(indel_dfs, folder_path + sample_ID + '_assembled_indel_dfs_' + date + '.pkl')

def main6(sample_ID, date, folder_path, cut_site_index):
	indel_dfs = retrieve_object(folder_path + sample_ID + '_assembled_indel_dfs_' + date + '.pkl')
	cut_loci = set(cut_site_index)
	temp = str(cut_site_index[0]) + '_' + str(cut_site_index[-1]) + '_'
	writer = pd.ExcelWriter(folder_path + sample_ID + '_assembled_' + temp + date + '.xlsx') 
	if 'wt_read' in indel_dfs.keys():
		df_wt = indel_dfs['wt_read']
		rev_index = []
		a = 0
		for i in df_wt.index:
			if int(i) > 200:
				rev_index.append(a)
			a += 1
		print rev_index
		df_wt_rev = df_wt.iloc[rev_index]
		df_wt_rev.to_excel(writer, 'wt_read')
	if 'deletion' in indel_dfs.keys():
		df_del = indel_dfs['deletion']
		rev_index = []
		a = 0
		for i in df_del.index:
			temp = set(list(n for n in range(i[0][0], i[0][1]+1)))
			if cut_loci & temp:
				rev_index.append(a)
			a += 1
		print rev_index
		df_del_rev = df_del.iloc[rev_index]
		df_del_rev.to_excel(writer, 'deletion')
	if 'insertion' in indel_dfs.keys():
		df_ins = indel_dfs['insertion']
		rev_index = []
		a = 0
		for i in df_ins.index:
			temp = set(list(n for n in range(i[0][0], i[0][1]+1)))
			if cut_loci & temp:
				rev_index.append(a)
			a += 1
		print rev_index
		df_ins_rev = df_ins.iloc[rev_index]
		df_ins_rev.to_excel(writer, 'insertion')
	if 'indel' in indel_dfs.keys():
		df_indel = indel_dfs['indel']
		rev_index = []
		a = 0
		for i in df_indel.index:
			temp = set(list(n for n in range(i[0][0], i[0][1]+1)))
			if cut_loci & temp:
				rev_index.append(a)
			a += 1
		print rev_index
		df_indel_rev = df_indel.iloc[rev_index]
		df_indel_rev.to_excel(writer, 'indel')					
	if 'mismatch' in indel_dfs.keys():
		df_indel = indel_dfs['mismatch']
		rev_index = []
		a = 0
		for i in df_indel.index:
			temp = set(list(n for n in range(i[0][0], i[0][1]+1)))
			if cut_loci & temp:
				rev_index.append(a)
			a += 1
		print rev_index
		df_indel_rev = df_indel.iloc[rev_index]
		df_indel_rev.to_excel(writer, 'mismatch')
	if 'multi_indel' in indel_dfs.keys():
		df_multi = indel_dfs['multi_indel']
		rev_index = []
		a = 0
		for i in df_multi.index:
			b = 0
			for j in i:
				for k in range(len(j)):
					temp = set(list(n for n in range(j[k][0][0], j[k][0][1]+1)))
					if cut_loci & temp:
						b += 1
			if b > 1:
				rev_index.append(a)
			a += 1
		print rev_index
		df_multi_rev = df_multi.iloc[rev_index]
		df_multi_rev.to_excel(writer, 'multi_indel')


if __name__ == '__main__':
	sample_ID = 'JG02-Fig1-N-2'
	date = '02012019'
	template_file = '/Users/leetz/Desktop/Jorge_papers_2017/CLADES_FigS1_12032018/Jorge_reference_seq_07162018.txt'
	cut_site_index = [193, 194, 195, 196, 197]
	folder_path = '/Users/leetz/Desktop/Jorge_papers_2017/CLADES_FigS1_12032018/'
	fp = open(template_file)
	seq_ori = fp.read()
	seq = seq_ori.upper()
	seq_list = list(seq)
	seq1 = ''.join(a for a in seq_list if a.isalpha())
	print len(seq1)	
	main1(sample_ID, date, folder_path)
	main2(sample_ID, date, folder_path)
	main3(sample_ID, date, folder_path, seq1)
	main4(sample_ID, date, folder_path, seq1)
	main5(sample_ID, date, folder_path, seq1)
	main6(sample_ID, date, folder_path, cut_site_index)
