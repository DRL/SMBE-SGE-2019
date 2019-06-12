#!/usr/bin/env python

import sys, itertools
import msprime
import numpy as np
import random


def main():
    index = sys.argv[1]

    directory = '/exports/csce/eddie/biology/groups/lohselab/sims/output/msprime/current_winner/'
    #directory = '/users/s1854903/sims/scripts/HeliconiusPopHist/msprime/test/'
    num_replicates = 1

    sample_size = 5

    mig = 3.8866e-7
    seqLength = 32e3 
    recr = 1.84675e-8
    Ne0 = 2.3241e6 #Cydno = ancestral
    Ne1 = 9.8922e5 #Mel_rosina = derived
    splitT = 4.8580e6
    secT = 1e5
    proportion = 0.1
    mu = 1.9e-9

    population_configurations = [
    msprime.PopulationConfiguration(sample_size=sample_size, initial_size=Ne0),
    msprime.PopulationConfiguration(sample_size=sample_size, initial_size=Ne1),
    ]
    
    #demographic events: specify in the order they occur backwards in time
    demographic_events = [
    msprime.MassMigration(time=secT, source=1, destination=0, proportion=proportion),
    msprime.PopulationParametersChange(time=splitT, initial_size=Ne0, population_id=0),
    msprime.MassMigration(time=splitT, source=1, destination=0, proportion=1.0)
    ]

    
    
    seed = np.random.randint(1, 2**32 - 1, num_replicates)

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
        with open(directory+'sim{}.vcf'.format(str(index)), 'w') as vcf_file:
            ts.write_vcf(vcf_file, ploidy=2)
        ts.dump(directory+'sim{}.trees'.format(str(index)))


if __name__ == '__main__':
    main()
