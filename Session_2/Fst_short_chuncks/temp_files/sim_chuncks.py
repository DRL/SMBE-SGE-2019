#!/usr/bin/env python

import sys, itertools, os
import msprime
import numpy as np
import random
from multiprocessing import Pool
from tqdm import tqdm

def run_simulation(params):
    
    seqLength, recr, population_configurations, demographic_events, mig, mu, seed, index = params
    #migration is asymmetric from pop_0 to pop_1, see migration_matrix
    replicates = msprime.simulate(
        num_replicates = 1,
        length = seqLength, 
        recombination_rate = recr,
        population_configurations = population_configurations,
        demographic_events = demographic_events,
        migration_matrix = [[0,0],
                            [mig,0]],
        mutation_rate = mu,
        random_seed=seed)
    for ts in replicates:
        msprime.mutate(ts, rate=mu, keep=True)
        with open('sim{}.vcf'.format(str(index)), 'w') as vcf_file:
            ts.write_vcf(vcf_file, ploidy=2)
        ts.dump('sim{}.trees'.format(str(index)))

    return True

def main():
    
    n = 30 #num_cores
    num_replicates = 30

    sample_size = 20

    mig = 3.8866e-7
    seqLength = 32e3 
    recr = 1.84675e-8
    Ne0 = 2.3241e6 
    Ne1 = 9.8922e5
    splitT = 4.8580e6
    mu = 1.9e-9

    population_configurations = [
    msprime.PopulationConfiguration(sample_size=sample_size, initial_size=Ne0),
    msprime.PopulationConfiguration(sample_size=sample_size, initial_size=Ne1),
    ]
    
    #demographic events: specify in the order they occur backwards in time
    demographic_events = [
    msprime.PopulationParametersChange(time=splitT, initial_size=Ne0, population_id=0),
    msprime.MassMigration(time=splitT, source=1, destination=0, proportion=1.0),
    msprime.MigrationRateChange(splitT, 0)
    ]
    
    params = [seqLength, recr, population_configurations, demographic_events, mig, mu]
    seeds = np.random.randint(1, 2**32 - 1, num_replicates)
    index_params = [params[:]+[seed, index] for index, seed in enumerate(seeds)]
    finished = []

    if n > 1:
        with Pool(n) as pool:
            with tqdm(index_params, total=len(index_params), desc="[%]", ncols=100) as pbar:
                for result in pool.imap_unordered(run_simulation, index_params):
                    finished.append(result)
                    pbar.update()    
    else:
        for prm in index_params:
            finished.append(run_simulation(prm)) 

if __name__ == '__main__':
    main()
