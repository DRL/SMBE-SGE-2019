import numpy as np
import pandas as pd
import msprime, tskit
import itertools, random, math
import threading
import tqdm
import intervals, intervaltree

#change weighting of p1,p2, between!!!
def calculate_Fst(ts, samplesize_p1=None, samplesize_p2=None, pop_id = (0,1), windowsize = 500e3, overlap=10e3, output='Fst0.csv'):
	#determine names of the different subsamples
	windowsize, overlap =int(windowsize), int(overlap)
	if not(samplesize_p1 and samplesize_p2):    
		samplesize_p1 = ts.get_sample_size()//2
		samplesize_p2 = ts.get_sample_size()//2
	else:
		subsample_p1 = random.sample(list(ts.get_samples(population_id=pop_id[0])), samplesize_p1)
		subsample_p2 = random.sample(list(ts.get_samples(population_id=pop_id[1])), samplesize_p2)
	
		subsample_combo = subsample_p1 + subsample_p2
		#treeSeq can be simplified to the tree for the selected samples
		ts = ts.simplify(subsample_combo)
		#resample tree, nodes have changed names    
	subsample_p1 = ts.get_samples(population_id=pop_id[0])
	subsample_p2 = ts.get_samples(population_id=pop_id[1])
	combinatorial_weights = (num_combos(samplesize_p1), num_combos(samplesize_p2), samplesize_p1*samplesize_p2)
	##########

	#calculating Fst
	#make dict containing list of windows and the relative weights of each tree across windows
	window_dict = windows(ts, windowsize, overlap)
	#calculate tree-wise Fst values
	tree_Fst_dict = tree_Fst(ts, subsample_p1, subsample_p2, combinatorial_weights)
	#combine tree-wise Fst values into window-wise values
	window_Fst_dict = window_Fst(tree_Fst_dict, window_dict)
	###########
	
	#generate output: dict to pandas dataframe
	return list(window_Fst_dict.values())

def windows(ts, windowsize, overlap):
	'''
	returns dict 
	key : window index
	value: tuple ([start, end window], [tree indexes of trees flanking window], [sequence length based weights
	of each of the trees])
	'''

	total = sum(1 for i in range(int(ts.get_sequence_length())//overlap) 
			if i*overlap+windowsize - 1 < ts.get_sequence_length())
	progress_bar = tqdm.tqdm(total=total)
	progress_bar.set_description(desc='determining windows')

	t_intervals = [tree.interval for tree in ts.trees()]
	tree_tree = intervaltree.IntervalTree.from_tuples(t_intervals)
	seq_len = int(ts.get_sequence_length())
	result_dict = dict()

	for index in range(seq_len//overlap):
		start, stop = index*overlap, index*overlap+windowsize-1  
		if stop < seq_len:
			temp_tree = tree_tree.copy()
			temp_tree.slice(stop)
			temp_tree.slice(start)
			window_breakpoints = sorted(temp_tree.overlap(start, stop))
			weights = [item.length()/windowsize for item in window_breakpoints]
			flanking = (ts.at(start).index, ts.at(stop).index)
			result_dict[index] = ((start, stop), flanking, weights)
			progress_bar.update()

	progress_bar.close()
	return result_dict

def num_combos(samplesize):
	return math.factorial(samplesize)/(math.factorial(samplesize-2)*2)

def tree_Fst(ts, subsample_p1, subsample_p2, combinatorial_weights):
	
	tree_Fst_dict = dict()
	num_threads = ts.get_num_trees()
	progress_bar = tqdm.tqdm(total=num_threads)
	progress_bar.set_description(desc='calculating tree-wise Fst')

	def thread_worker(thread_index):
		#for each tree calculate Fst
		tree = ts.at_index(thread_index)
		p1 = np.mean([tree.tmrca(u,v) for u,v in itertools.combinations(subsample_p1, 2)])
		p2 = np.mean([tree.tmrca(u,v) for u,v in itertools.combinations(subsample_p2, 2)])
		within_av = np.mean([p1,p2])
		between = np.mean([tree.tmrca(u,v) for u,v in itertools.product(subsample_p1, subsample_p2)])
		data_av = np.mean([within_av,between])
		tree_Fst_dict[thread_index] = (data_av - within_av)/data_av
		progress_bar.update()

	threads = [
		threading.Thread(target=thread_worker, args=(i,)) 
			for i in range(num_threads)]
	for t in threads:
		t.start()
	for t in threads:
		t.join()
	progress_bar.close()
	return tree_Fst_dict

def window_Fst(tree_Fst_dict, window_dict):
	window_Fst_dict = dict()
	num_windows = len(window_dict)
	progress_bar = tqdm.tqdm(total=num_windows)
	progress_bar.set_description(desc='recombining tree-wise into window-wise Fst values')

	def thread_worker(thread_index):
		#for each window recombine tree-wise Fst values into window-wise values
		window, indices, weights = window_dict[thread_index]
		start, stop = window
		window_Fst_dict[thread_index] = (start, stop, np.average([tree_Fst_dict[index] 
			for index in range(indices[0], indices[1]+1)], weights=weights, axis=0))
		progress_bar.update()

	threads = [
		threading.Thread(target=thread_worker, args=(i,)) 
			for i in range(num_windows)]
	for t in threads:
		t.start()
	for t in threads:
		t.join()
	progress_bar.close()

	return window_Fst_dict