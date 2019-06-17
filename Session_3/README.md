# Session 3

# ```gIMble``` part 1: Filtering, blocking and visualising

## 0. Introduction

The software ```gIMble``` (**g**enome-wide **IM** **b**lockwise **l**ikelihood **e**stimation) is a toolkit for compressing 
coverage data – in the form of BED intervals shared across samples – and
variation data – as records in a VCF file – into blocks of uniform length and constant sample coverage for which we record certain types of mutation configurations. 

These blocks are the basic unit of data for subsequent demographic inference.     

### 0.1 Input data
We will analyse variation within and between two species of *Heliconius* butterflies – *H. melpomene rosina* (‘ros’) and *H. cydno chioneus* ('chi') – along chromosome 18 ('Chr18') of the reference genome ```Hmel2``` (Davey et al, 2017). Each population is composed of 10 samples. All input files for this session can be found in the folder ```/data/hmel.chr18/```. We will use the input data (VCF file, BED multinter coverage file, genome file and sample file) to construct 'blocks' of 64 b across the genome, record the variation within the blocks across the samples that cover it, filter the blocks and port them to an alternative coordinate system, and group blocks into sliding windows across the genome to compute genome scans and generate the input data for the demograghy analysis we will tackle in the following session.      

### 0.2 Input data formats

**Variation data**
- Variation data is often encoded in [VCF](http://www.internationalgenome.org/wiki/Analysis/vcf4.0/) (Variant Call Format) files. 
- VCF files are the result of DNA sequencing data being subjected to a complex analytical (and often, non-linear!) process composed of read mapping against a reference genome, variant calling and data QC. 
- These files are usually compressed using ```bgzip``` and indexed via ```tabix``` to ease handling and access. For more information, visit http://www.htslib.org/.   
- Typically, only sites that show variation in at least one sample are recorded. This is based on the assumption that all sites in the reference were uniformly covered in all sequenced samples. In reality, sequencing coverage of genomic regions tends to vary stochastically between samples ![Figure 1](https://github.com/DRL/SMBE-SGE-2019/blob/master/Session_3/figure_1.png "Figure 1")


This has obvious implications for population genetic inference since it blurs the distinction between variable and monomorphic sites.

**Coverage data**
- One way of dealing with the issue of differential coverage is identifying genomic regions with acceptable coverage (above/below certain thresholds) in each sample. This can be achieved using [GATK’s CallableLoci Tool](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_coverage_CallableLoci.php) – to identify “callable” regions based on the BAM files for each sample – followed by intersection of the resulting BED files using [bedtools multiinter](https://bedtools.readthedocs.io/en/latest/index.html). This yields a BED file in which the last column of each interval describes the set of samples that cover it at acceptable coverage, e.g.:

```
chrom	start	end	num	list
chr1	0	100	6	A,B,C,X,Y,Z
chr1	100	164	4	A,B,X,Y
chr1	164	228	4	A,C
chr1	228	267	4	Y,Z
chr1	292	356	4	B,X,Y,Z
chr1	356	410	4	Z
[...]
```

**Genomic sequence data**
- While the BED file lists callable regions across the genome, it lacks information about sites not covered in any sample. This is a problem once one wants to plot metrics along the genome, since the software has to know how long each sequence is and where the data points need to be positioned. For this purpose we are using genome files as used in BedTools composed of two columns containing sequence ID and length (see [genome files](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#genome-file-format)), e.g.:

```
chr1	305046
chr2	47928
chr3	1895173
[...]
```

**Sample information data**
Analysis of variation within and between samples of two populations/species is only possible if the samples are labelled by population. We do this using a sample file (CSV) in the following format:

```
A,pop1
B,pop1
C,pop1
X,pop2
Y,pop2
Z,pop2
```

### 0.3 Data consistency
Successful execution of ```gIMble``` analyses is only possible if variation, coverage, genomic sequence and sample information data can be linked within gIMble. Therefore, the SAME sequence/sample IDs have to be used across the files.

### 0.4 Output data formats
Output files written by ```gIMble``` are either PNG (Portable Network Graphics) image files or tables in HDF (Hierarchical Data Format) files. Tables from the HDF output files can easily be exported into CSV/TSV files for post-processing by using the script ```hdf5v.py``` in the ```~/gIMble``` folder.

## 1. ```gIMble``` modules
In this session, we are going to analyse the ```Heliconius``` data for Chr18 using five ```gIMble``` modules:

* ```gIMble blocks```: constructs blocks (of fixed length) along the genome, while maximising coverage across pairs of samples between the two species based on the BED file.

* ```gIMble variants```: queries genotypes from the VCF file for the samples present in each block and records block-wise variation data. 

* ```gIMble modify blocks```: transforms the coordinate system from *Heliconius* ```Hmel``` scaffolds/contigs to ```Hmel``` chromosomes using a coordinate mapping file.

* ```gIMble modify variants```: filters variation data based on the occurrence of missing/multiallelic genotypes within blocks. 

* ```gIMble windows```: constructs sliding windows of constant number of blocks across the genome and records window-wise variation data.

### 1.0 Concepts

In order to understand what is happening under the hood when running ```gIMble``` it is necessary to introduce the following concepts:

* **BED interval**: a continuous region in the BED file for which a set of samples display adequate (i.e. 'callable') coverage.
 
* **Interspecies pair (IP)**: all combinations of sample IDs between the two species specified in the sample file given to ```gIMble```. Based on these, BED intervals are being considered for inclusion in blocks. The ```Heliconius``` sample file includes 10 samples for each species, hence a given block can have up to 100 possible IPs.

* **block**: basic unit of analysis used within ```gIMble```. A block has the following properties:
  * It is composed of BED intervals 
  * Its length – the sum of lengths of its BED intervals – is uniform and specified at the beginning of the analysis. Based on the uniform length we can assume &theta; to be comparable across blocks. (During this practical we will construct blocks that are 64 b long).
  * The IP composition of a block is the intersection of IPs across its BED intervals. For all its IPs – and the underlying samples – both monomorphic and variant sites are recorded.
  * A block can contain gaps, as long as its span – the distance between the start of its first and the end of its last BED interval – does not exceed a given threshold. We will set this threshold to 80 b, thereby allowing up to 16 b between BED intervals to not be covered by any of the block’s IPs. 
  * Variation within a block is recorded for each IP and is encoded as a **mutuple** consisting of four categories (mutypes): **hetA**, **fixed**, **hetB**, and **hetAB**. For some examples as to how genotypes are encoded as mutuples see the following table. In addtion, two other counts per block and IP are recorded – the number of multiallelic and missing genotypes; based on these blocks can be included/excluded from the analysis.

| GT of sample A  | GT of sample B |Mutype| 
| ------------- | ------------- | ------------- |
| 0/1  | 0/0  |	hetA |
| 1/1  | 0/0  |	fixed |
| 0/0  | 2/2  |	fixed |
| 0/0  | 0/1  |	hetB |
| 0/1  | 0/1  |	hetAB |
| 1/2  | 1/2  |	hetAB |
| 0/1  | 1/2  |	multiallelic (non-biallelic site) |
| -/-  | 0/1  |	missing |

* **window**: a window is composed of a set of neighbouring blocks with uniform count. Variation within windows is recorded as the tally of mutuples contained in its blocks.

## 1.1 ```gIMble blocks```

We will construct blocks based on the BED intervals contained in the file ```hmel.chr18.bed```. We will pass sample and genome sequence length information via the files ```hmel.samples.csv``` and ```hmel.chr18.genomefile```, respectively. We also will run the analysis using four threads (```-t 4```).

In order to test the effect of coverage on the final results of our analyses we will run two independent analyses using different values for the minimum number of IPs that have to be present in a BED interval for it to be considered for block construction:

- ```-n 5``` : at least five samples in both populations have to be present in a BED interval  
- ```-n 10```: at least ten (i.e. all) samples in both populations have to be present in a BED interval

To be able to distinguish the output files, we will give them unique names by specifying output prefixes using the ```-o``` argument.

``` 
# n = 5
./gIMble blocks \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-b /data/hmel.chr18/hmel.chr18.bed \
	-o hmel.chr18.n_5 \
	-n 5 \
	-t 4
	
# n = 10
./gIMble blocks \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-b /data/hmel.chr18/hmel.chr18.bed \
	-o hmel.chr18.n_10 \
	-n 10 \
	-t 4
```

During the analysis, ```gIMble blocks``` will print status messages. Based on these, please try to answer the following questions:

**Questions**:
- In which way do the two analyses differ? 
- Which analysis was more successful? Why? Define 'successful'.
 
Three output files are generated by ```gIMble blocks```:
- ```*.blocks.h5``` : an HDF5 datastore containing information about the generated blocks. This file will be used in subsequent analyses steps.
- ```*.distance.png``` : a log-log scatter plot of the counts of distances between blocks. This information is relevant for subsequent processing of the blocks using ```gIMble windows``` (how far apart should blocks be to still end up in the same window?).
```*.blocks_per_sample.png``` : barcharts of number of bases in blocks by sample ID.

**Questions**:
- How do the ```*.distance.png``` plots differ?
- How do the ```*.blocks_per_sample.png``` plots differ?

Using the program ```hdf5v.py``` (located in the ```gIMble/``` folder), one can investigate the content of a HDF5 datastore, e.g.:

```
# display the names of tables contained in the datastore
./hdf5v.py -f hmel.chr18.n_5.blocks.h5 

# export the table ‘block’ as a TSV file (this will print to screen basic metrics of the numeric columns) 
./hdf5v.py -f hmel.chr18.n_5.blocks.h5  -d block -t 

# export all tables as CSV files
./hdf5v.py -f hmel.chr18.n_5.blocks.h5  -c 
```

## 1.2 ```gIMble variants```

For the blocks we created in the previous step, we will now query the VCF file for variation of the samples contained in the IPs. For this we have to pass the VCF file ```hmel.chr18.vcf.gz``` using the ```-v``` argument and the previously created blocks HDF5 datastore.

Again, we will pass sample and genome sequence length information via the files ```hmel.samples.csv``` and ```hmel.chr18.genomefile```, respectively. We also will run the analysis using four threads (```-t 4```). Since this step is quite time intensive, it is probably best to run only one of the two analyses (n=5 **or** n=10). Try coordinating with the people around you to divide the tasks. While this analysis is running it might also be a good opportunity to ask any questions you might have.

```
# n = 5
./gIMble variants \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-v /data/hmel.chr18/hmel.chr18.vcf.gz \
	-o hmel.chr18.n_5 \
	-b hmel.chr18.n_5.blocks.h5 \
	-t 4

# n = 10
./gIMble variants \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-v /data/hmel.chr18/hmel.chr18.vcf.gz \
	-o hmel.chr18.n_5 \
	-b hmel.chr18.n_5.blocks.h5 \
	-t 4
```

The module ```gIMble variant``` will print status messages to screen as well as metrics for samples, populations, and the entire dataset. 

**Questions**:
- Do the two analyses differ in their amount of monomorphic blocks? Why is that?

The sample-wise metrics describe the number of sites in blocks for each of the samples, as well as their heterozygosity (count of heterozygous sites divided by the number of total sites) and proportion of sites with missing genotypes in blocks. The population-wise metrics contain the averages of these values across their respective samples. The dataset metrics, however, are calculated based on the mutuples (mutation configurations) in all blocks across their IPs. FGVs (four-gamete-violations) is the proportion of blocks displaying both ‘fixed’ and ‘hetAB’ mutypes, which suggest incompatible underlying trees. These blocks are not used when computing metrics and are a proxy for dataset quality. High FGV proportions in a dataset suggest that the selected block length is too long. 

Metrics for &pi;, *D*<sub>xy</sub>; and *F*<sub>st</sub>; are based on all non-FGV blocks in the dataset. The reason for &pi; of a population being lower than its heterozygosity, resides in the fact that mutuples of well covered, conserved (i.e. monomorphic) blocks contribute more than those of poorly covered ones.

Two output files are generated:
- ```*.variants.h5```: the HDF5 datastore of variation data of blocks
- ```*.mutuple_barchart.png```: barcharts of the counts/proportion of the most frequent mutypes in the dataset.

## 1.3 ```gIMble modify```  
In this and the next step we will modify the existing block and variant HDF5 datastores using the ```gIMble modify``` module. With ```gIMble modify blocks``` we can transfer the blocks from one coordinate system to another and with ```gIMble modify variants``` we can exclude mutuples from blocks based on their count of missing/multiallelic genotypes.  

```
# n = 5
./gIMble modify blocks \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-c /data/hmel.chr18/hmel.chr18.chrom_coordinates.txt \
	-o hmel.chr18.n_5 \
	-b hmel.chr18.n_5.blocks.h5

# n = 10
./gIMble modify blocks \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-c /data/hmel.chr18/hmel.chr18.chrom_coordinates.txt \
	-o hmel.chr18.n_10 \
	-b hmel.chr18.n_10.blocks.h5
```

Each of these commands will generate three files (analogous to those of ```gIMble blocks```)
- ```*.modified.blocks.h5```
- ```*.modified.distance.png```
- ```*.modified.blocks_per_sample.png```

The plots should show no differences to the previous plots generated with ```gIMble blocks```, as all blocks are being retained and converted into the coordinate system specified in the file ```hmel.chr18.chrom_coordinates.txt``` (you can use ```less``` to actually look into the file). If some contigs in the dataset were not allocated chr18, these would be excluded from the output.

To filter mutuples associated with high counts of missing genotypes in blocks (we will not filter on multiallelic genotype counts), run one of the following commands based on which variant HDF5 datastore you have created previously.
 
```
# n = 5
./gIMble modify variants \
	-v hmel.chr18.n_5.variants.h5 \
	-o hmel.chr18.n_5 \
	-m 4 \
	-M 64

# n = 10
./gIMble modify variants \
	-v hmel.chr18.n_10.variants.h5 \
	-o hmel.chr18.n_10 \
	-m 4 \
	-M 64
```

This will generate output files analogous to the ```gIMble variant``` module. However, no metrics concerning samples, populations or the entire datasets are printed to screen.

**Questions**:
- How many mutuples are excluded based on the applied filter?
- Which mutuples are excluded? How can you find out?

## 1.4 ```gIMble windows```

In the last step we will join neighbouring blocks into sliding windows with a defined number of blocks and a defined overlap of blocks between windows. This will allow us to plot metrics across the genome and collect mutuple tallies for the windows which we will use in the inference session. Keep in mind that since we have changed the coordinate system we now have to provide the genomefile of chr18 ```hmel.chr18.new_coordinates.genomefile```. 

```
# n = 5
./gIMble windows \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.new_coordinates.genomefile \
	-b hmel.chr18.n_10.modified.blocks.h5 \
	-v hmel.chr18.n_10.modified.variants.h5 \
	-o hmel.chr18.n_10 \
	-w 500 \
	-l 100 \
	-t 4

# n -10
./gIMble windows \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.new_coordinates.genomefile \
	-b hmel.chr18.n_10.modified.blocks.h5 \
	-v hmel.chr18.n_10.modified.variants.h5 \
	-o hmel.chr18.n_10 \
	-w 500 \
	-l 100 \
	-t 4
```

This will create four output files:
- ```*.fst_genome_scan.png``` : a *F*<sub>st</sub>; scan across chr18
- ```*.pi_genome_scan.png``` : a &pi; scan across chr18
- ```*.pi_scatter.png``` : &pi; vs. *D*<sub>xy</sub>; scatter plots of all windows in chr18
- ```*.windows.h5``` : the HDF5 datastore of windows

**Questions**:
- How do the genome scans differ between the n=5 and n=10 datasets?
- Look again at the ```*.distance.png``` file you generated earlier. Could you use this plot to select a more appropriate value for the argument ```--max_block_distance``` in ```gIMble windows```? Do the results change?

In the next session we will use a global ```*.variants.h5``` file (computed on the entire ```Heliconius``` genome) to estimate parameters under the divergence model (using the ```gIMble likelihood``` module) and the ```*.windows.h5``` file to perform a gridsearch across the windows of chr18 (using the ```gIMble gridsearch``` module). 
