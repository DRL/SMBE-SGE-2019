# gIMble part 2: Fitting background demography

The aim of this practical is to infer demographic parameters under a simple set of demographic models. We will initially focus on a genome-wide background history and use a pre-computed global ```*.variants.h5``` file (computed on the entire ```Heliconius``` genome) containing the tallies of mutuples for the two *Heliconius* butterflies – *H. melpomene rosina* (‘ros’) and *H. cydno chioneus* ('chi'). This is analogous to the file we generated in the last session for chromosome 18 but for the whole genome.

## Input files
The data for which we want to estimate demographic parameters are the pre-computed global mutuple tally – the counts of distinct mutuples across the whole *Heliconius* genome (not including the Z chromosome) for the 20 samples in question – and is stored in the file ```~/gIMble/smbe/hmel.all_mutype_counts.h5```.

To be able to fit a demographic model we require a ‘model’ file, in this case```~/gIMble/models/divergence.txt```. In addition, we need to specify the sample (and genome) file we have used in the previous session. 

## ```gIMble likelihood```
The divergence model which we are using is simple (two parameters &theta; – assumed to be the same for both species and the ancestor for now – and *T*), i.e. one population is constrained to have the same *N*<sub>e</sub> as the ancestral population. 

![Divergence model](https://github.com/DRL/SMBE-SGE-2019/blob/master/Session_4/divergence.jpg "Divergence model")

```gIMble likelihood``` employs a heuristic simplex algorithm (Nelder-Mead) to maximise composite likelihood (lnCL). For the likelihood search, one can either fix certain parameters or provide upper and lower bounds which are then used in the starting simplex of the heuristic search. 

We will try to fit a simple divergence model by forcing *N*<sub>e</sub> to be the same for both ancestral and derived populations. We do this by labelling one of the populations as ancestral (```-A “ros”```) and setting *N*<sub>e</sub> of the derived population to 1.0 (```--derived_Ne 1.0```). 

To be able to deal with the data efficiently, we will only calculate probabilities of mutuples whose mutype counts are at most 2 (```-k 2```). Mutuples containing counts above 2 will be "binned" and their probabilities of occurrence on the state graph are calculated as marginal probabilities (e.g. [3, 0, 0, 0], [4, 0, 0, 0], [5, 0, 0, 0] … are treated as [≥3, 0, 0, 0]).
   
Based on the model, we will estimate the two parameters: &theta; (per block) and the divergence time *T* for which we set lower and upper boundaries:

```
./gIMble likelihood \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-l models/model.divergence.txt \
	-A "ros" \
	-v smbe/hmel.all_mutype_counts.h5 \
	-k 2 \
	--theta_low 0.5 \
	--theta_high 1.4 \
	--time_low 0.2 \
	--time_high 0.6 \
	--derived_Ne 1.0 \
	-e 12345 \
	-t 4 
```

## Model parameters estimated via Mathematica

| model | lnCL | &theta; | *C*<sub>derived</sub> | *T*|
| -----| ---| ----| ----| ----|
| Divergence | -4.9847270E+08 | 1.12614 | 1 | 0.3820622|
| Div2 (Ancestor = 'ros') | -4.9832980E+08 | 1.090111 | 0.8241344 | 0.4140925|
| Div2 (Ancestor = 'chi') | -4.9799600E+08 | 1.196783 | 1.401532 | 0.3229429|

## More complex models

We now want to allow for migration at a constant rate *M*=4\**N*<sub>e</sub>*m*. 

![IM model](https://github.com/DRL/SMBE-SGE-2019/blob/master/Session_4/IM.pdf "IM model")

However, maximising lnCL for the IM model takes several hours (and is not implemented stabily in gIMble likelihood) so we will simply examine precooked point estimates:

|Model | lnCL | &theta; | *C*<sub>derived</sub> | *T* | Migration|
| -----| ---| ----| ----| ----| --- |
|IM2 (Ancestor = 'ros', 'chi'->'ros') | -4.9833E+08 | 1.0100 | 0.7504 | 0.5516 | 0.7490|
|IM2 (Ancestor = 'chi', 'chi'->'ros') | -4.9700E+08 | 1.0873 | 2.5783 | 1.4999 | 4.0728|


|Model | lnCL |&theta; | *C*<sub>derived</sub> | *T* | Migration |
| -----| ---| ----| ----| ----| --- |
|IM2 (Ancestor = 'ros', 'ros'->'chi') | -4.9822E+08 | 1.0643 | 1.0408 | 0.5323 | 1.2735|
|IM2 (Ancestor = 'chi', 'ros'->'chi') | -4.9799E+08 | 1.1904 | 1.3832 | 0.3351 | 0.1920|

**Questions**: 
* There are four possible IM models: Which one fits best and how do *T* and *N*<sub>e</sub> differ from the divergence models we fitted previously? 
* By how much does model support increase when we include gene flow? 
* What is the probability of migration per lineage and generation (this is what we would need to set up an msprime simulation)?

## The effect of monomorphic blocks

Monomorphic blocks (0,0,0,0) are the most frequent class of blocks and, unsurprisingly, have a disproportionate effect on the inference. We can remove this effect by conditioning on only observing blocks that contain variants (still to be implemented in gIMble).

### Divergence model
| Model | lnCL | &theta; | *C*<sub>derived</sub> | *T* |
| -----| ---| ----| ----| ----| 
| Div1 | -4.310784E+08 | 1.24152 | 1.00000 | 0.39105|
| Div2 (Ancestor = 'ros') | -4.307465E+08 | 1.18526 | 0.74536 | 0.44068|
| Div2 (Ancestor = 'chi') | -4.308385E+08 | 1.29252 | 1.26807 | 0.34838|

### IM model
|Model | lnCL | &theta; | *C*<sub>derived</sub> | *T* | Migration|
| -----| ---| ----| ----| ----| --- |
|IM2 (Ancestor = 'ros', 'chi'->'ros') | -4.2965E+08 | 0.8862 | 0.5387 | 1.0780 | 1.4282|
|IM2 (Ancestor = 'chi', 'chi'->'ros') | -4.2988E+08 | 1.1990 | 2.1108 | 0.8191 | 3.1439|

# gIMble part 3: Slippin' and slidin'

How do migration and *N*<sub>e</sub> vary along the genome? We could try and fit the IM model for each window along the genome. However, given the search time, this is clearly inefficient. A better alternative is to precompute tables of blockwise probabilities over a grid in parameter space. While a ```gIMble gridsearch``` would be inefficient for a single optimisation it becomes very efficient if we want to optimise the same model for many different datasets, in this case counts of blockwise configs for many sliding windows. Other Inferences such as ```sweepfinder``` and ```volcanofinder``` also rely on pre-computed grids. 

We will restrict ourselves to chromosome 18 and use a precomputed grid (18x12X12) for the best fitting IM2 model (Ancestor = 'ros', 'chi'->'ros') for migration, &theta;<sub>Ancestor</sub> and &theta;<sub>Derived</sub>, while fixing *T* to its global mCLE (Maximum composite Likelihood Estimate) .  

## Input parameters

We want to estimate demographic parameters along Chr18, hence the data on which we will run the ```gIMble gridsearch``` is the window-wise mutuple tallies (```hmel.chr18.n_5.windows.h5```) we have generated earlier. The model we will specify is located in ```models/model.IM.M_D2A.MM_D2A.txt``` and the grid of probabilities of mutuples by parameter set can be found in ```grid/grid.csv```. Furthermore, we need to specify the mutation rate (1.9e-9) as well as the length of blocks used in the analysis to be able to scale parameters as absolute estimates. Split time *T* is fixed to its global mCLE of 4e6.

## ```gIMble gridsearch```
```
./gIMble gridsearch \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.new_coordinates.genomefile \
	-l models/model.IM.M_D2A.MM_D2A.txt \
	-A "ros" \
	-w hmel.chr18.n_5.windows.h5 \
	-k 2 \
	-o hmel.chr18.min_10_sample \
	--mu 1.9e-9 \
	--block_size 64 \
	--time_MLE 4e6 \
	--grid grid/grid.csv \
	-t 4 
```

This will generate two output files:
* ```*.composite_likelihoods.h5```: HDF5 datastore containing parameter estimates
* ```*.parameter_scan.png```: a genome scan plot of &theta;<sub>Ancestor</sub>, &theta;<sub>Derived</sub> and migration across Chr18

We are going to release a new ```gIMble``` version soon which will allow construction of custom, user-specified grids to simplify parameter scans across genomes.

# Discussion questions:
* Why should *T* be fixed globally? 
* Are there alternative cartoon models that would be more useful/insightful for a sliding window analysis?
* How could we incorporate heterogeneity in recombination rate in this inference?
* Is there a simulation based check for this inference?
