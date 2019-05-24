# Session 3

# 0. Preparation 

- Conda environment: 
  - gimble.sage.yaml (hopefully not needed)

- Data:
  - precomputed Heliconius data in /data/hmel
  - example dataset in /data/example

- Other
  - add gIMble executable to path 

# 1. Generate state graph for different models \[optional\] (TIMING)
```
./gIMble graph -s input/test.minimal.samples.shuffled.csv [MODELFLAGS] -o output/gimble
```

## Inspect results
- State graphs (pngs)
- Model files (tsv)

# 3. Calculate PODs with fixed parameters (TIMING)
```
python2 src/probs.py pods -p output/master.2_pop.2_ploidy.E.paths.txt -t 1 -A 1.0 -D 1.0 -M 2.34 -m 1.2 -T 1.4 
```

## Inspect results
- POD file

# 3. Estimate parameters based on global bSFS (TIMING)
```
./gIMble estimate -b output/gimble.global_bSFS.txt -s input/test.minimal.samples.shuffled.csv -o output/gimble [BOUNDARIES]
```

# Inspect results
- Summary file (text)

# 4. Estimate parameters based on window bSFS (TIMING)
```
./gIMble estimate -b output/gimble.window_bSFS.txt -s input/test.minimal.samples.shuffled.csv -o output/gimble [BOUNDARIES]
```

# Inspect results
- Summary file (text)
- parameters by windows (png)

# 5. Inspect Heliconius data results
...