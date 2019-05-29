import msprime, tskit
import numpy as np

def determine_windows(ts, windowsize, overlap):
	if overlap:
		assert windowsize%overlap==0, 'windowsize is not a multiple of the overlap'
		windows = [i*overlap for i in range(int(ts.sequence_length)//overlap)+1]
	else:
		windows = [i*windowsize for i in range(int(ts.sequence_length)//windowsize)+1]
	
	if windows[-1] != ts.sequence_length:
			windows.append(int(ts.sequence_length))

def Fst(div_entry):
	A = np.mean((div_entry.diagonal()))
	B = div_entry[1,0]
	return (B-A)/(A+B)

def calculate_Fst(ts, samplesets, windowsize, overlap=None):
	'''
	samplesets = [p1, p2]

	'''
	windows = determine_windows(ts, windowsize, overlap)
	div_matrix = branchStat.divergence_matrix(samplesets, windows)
	if overlap:
		jump = windowsize//overlap + 1
		div_matrix = np.mean([div_matrix[i:i+6] for i in range(len(windows)-6)],axis=1)

	Fst = [Fst(entry) for entry in with_overlap_mean]

	return Fst