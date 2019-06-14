import msprime, tskit
import numpy as np
import pandas as pd

def determine_windows(ts, windowsize, overlap):
	
	assert windowsize <= ts.sequence_length, 'windowsize is too large for your treesequence'
	if windowsize == ts.sequence_length:
		return [0, ts.sequence_length]
	else:
		if overlap:
			assert windowsize%overlap==0, 'windowsize is not a multiple of the overlap'
			windows = [i*overlap for i in range(int(ts.sequence_length)//overlap +1)]
		else:
			windows = [i*windowsize for i in range(int(ts.sequence_length)//windowsize + 1)]
		
		if windows[-1] != ts.sequence_length:
				windows.append(int(ts.sequence_length))
		return windows

def Fst(div_entry):
	A = np.mean((div_entry.diagonal()))
	B = div_entry[1,0]
	return (B-A)/(A+B)

def alt_Fst(div_entry):
	A = np.mean((div_entry.diagonal()))
	B = div_entry[1,0]
	return 1- (A/B)

def calculate_Fst(ts, samplesets, windowsize, overlap=None, function=Fst):
	'''
	samplesets = [p1, p2]

	'''
	windows = determine_windows(ts, windowsize, overlap)
	BranchStat = tskit.BranchLengthStatCalculator(ts)
	div_matrix = BranchStat.divergence_matrix(samplesets, windows)
	
	if not overlap or overlap==0:
		start, stop = [[0], [ts.sequence_length]]

	else:
		jump = windowsize//overlap + 1
		div_matrix = np.mean([div_matrix[i:i+jump] for i in range(len(windows)-jump)],axis=1)
		start = [windows[i] for i in range(len(windows)-(jump+1))]
		stop = [windows[i]+windowsize for i in range(len(windows)-(jump+1))]

	Fst_result = list(map(function, div_matrix))
	data_tuples = list(zip(start, stop, Fst_result))
	data = pd.DataFrame(data_tuples, columns=['start', 'stop', 'Fst'])
	
	return data