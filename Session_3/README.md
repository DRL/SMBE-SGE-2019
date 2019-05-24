# Session 3

# 0. Preparation 

- Conda environment: 
  - gimble.env.yaml (subject to change)
  - gimble.sage.yaml (hopefully not needed)

- Data:
  - precomputed Heliconius data in /data/hmel
  - example dataset in /data/example

- Other
  - add gIMble executable to path 

# 1. Generate blocks (TIMING)
```
./gIMble blocks -s input/test.minimal.samples.shuffled.csv -b input/test.minimal.bed -g input/test.minimal.genomefile -o output/gimble
```

# Inspect results
- Summary file (text)
- Counts of blocks per sample (png)
- Counts of blocks per sequence (png)
- Percentage of genome in blocks by pairs (png) 

# 2. Fetch variants (TIMING)
```
./gIMble variants -s input/test.minimal.samples.shuffled.csv -v input/test.minimal.vcf.gz -b output/gimble.blocks.bed -o output/gimble
```

# Inspect results
- Summary file (text)
	- monomorphic blocks
	- dxy/fst/pi per sample/pair
- Missing genotypes in Blocks (png)
- Multiallelic genotypes in Blocks (png)

# 3. Filter/fixcoordinates (TIMING)
```
./gIMble filter -b output/gimble.blocks.bed -s input/test.minimal.samples.shuffled.csv -c input/coordinates_map.tsv -o output/gimble
```

# Inspect results
- Summary file (text)
- Counts of blocks per sample (png)
- Counts of blocks per sequence (png)
- Percentage of genome in blocks by pairs (png) 

# 4. Windows (TIMING)
```
./gIMble windows -b output/gimble.blocks.bed -s input/test.minimal.samples.shuffled.csv -c input/coordinates_map.tsv -o output/gimble
```

# Inspect results
- Summary file (text)
- Counts of blocks per windows (png)
- Counts of windows per sequence (png)
- dxy/fst/pi of windows (png)

# 5. Inspect Heliconius data results
...