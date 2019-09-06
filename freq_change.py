import fwdpy11
import sys
import gzip
import pandas as pd
import sqlite3
from collections import namedtuple
import numpy as np

AlleleFreq = namedtuple(
    'AlleleFreq', ['generation', 'position', 'origin', 'dac',
                   'fixation', 'neutral', 'label', 's', 'w'])

PopFitness = namedtuple(
    'PopFitness', ['generation', 'mean', 'variance'])

def count_frequencies(pop, num_neutral):
    N = pop.N
    fixations = []
    for i, j in zip(pop.fixation_times, pop.fixations):
        if i >= 10*pop.N and j.g <= 10*pop.N + 200:
            fixations.append((j.g, j.pos))

    if num_neutral > 0:
        # Add neutral fixations, if there are any
        # NOTE: we have no idea about fixation times
        # for neutral mutations, as they were not simulated!
        for i, j in enumerate(pop.mcounts):
            if j == 2*pop.N and pop.mutations[i].neutral is True:
                # Neutral mutations have an origin time
                # randomly assigned based on branch lengths.
                # It is just uniform along the branch.
                g = pop.mutations[i].g
                if g >= 10*pop.N and g <= 10*pop.N + 200:
                    fixations.append((g, pop.mutations[i].pos))

    data = []
    pop_data = []
    # Iterate over all samaple nodes.
    # False means to exclude the "alive" individuals
    # at time pop.generation
    for t, n, metadata in pop.sample_timepoints(False):
        if t%N == 0:
            print(t)
        tables, idmap = fwdpy11.simplify_tables(pop.tables, n)
        muts_visited = 0

        # Convert input (pre-simplification) node labels
        # to output (post-simplification) node labels
        output_nodes = idmap[n]
        
        # we remap the input/output nodes to get the leaves at a given timepoint that carries a mutation
        # below in ws_samples
        remap_nodes = []
        for i in range(len(metadata)):
            remap_nodes.append(idmap[metadata[i]['nodes']])
        
        remap_nodes = np.ravel(remap_nodes)
        
        # Get the distribution of diploid fitnesses for this generation
        ws = metadata['w']
        pop_data.append(PopFitness(t, ws.mean(), ws.var()))
        
        # NOTE: do not update samples lists below nodes.
        for tree in fwdpy11.TreeIterator(tables, idmap[n], update_samples=True):
            for m in tree.mutations():
                assert m.key < len(pop.mutations)
                muts_visited += 1
                mut_node = m.node
                # NOTE: infinite sites means
                # leaf counts are the frequencies
                dac = tree.leaf_counts(mut_node)
                pos = tables.sites[m.site].position
                assert pos == pop.mutations[m.key].pos
                assert pos >= tree.left and pos < tree.right
                origin = pop.mutations[m.key].g
                fixed = int((origin, pos) in fixations)
                
                # Get the ws for the samples carrying the mutation
                m_samples = tree.samples_below(m.node)
                for mi in m_samples:
                    assert mi in output_nodes, f"{mi} not an output node"
                # There are two nodes per individual, hence the //2
                ws_samples = [metadata[np.where(remap_nodes == samp)[0][0]//2]['w'] for samp in m_samples]
                w_m = np.mean(ws_samples)
                
                # Edge case: we skip mutations that fixed prior
                # to any of these time points
                fixed_before = False
                if dac == 2*pop.N and fixed == 0:
                    fixed_before = True
                if not fixed_before:
                    data.append(AlleleFreq(t, pos, origin,
                                           dac, fixed, int(m.neutral),
                                           pop.mutations[m.key].label,
                                           pop.mutations[m.key].s, w_m))
        assert muts_visited == len(tables.mutations)
    return data, pop_data


if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]
    outfile_pop = sys.argv[3]
    neutral_mu = float(sys.argv[4])
    seed = int(sys.argv[5])

    pop = fwdpy11.DiploidPopulation.load_from_file(infile)

    rng = fwdpy11.GSLrng(seed)
    print(len(pop.mcounts))
    # NOTE: this function initiates mutation counting
    # via tree sequences as a side-effect. :(
    n = fwdpy11.infinite_sites(rng, pop, neutral_mu)
    print(len(pop.mcounts))
    print(f"{n} neutral variants added")
    
    #freqs = count_frequencies(pop, n)
    freqs, pop_data = count_frequencies(pop, n)
    df = pd.DataFrame(freqs, columns=AlleleFreq._fields)
    df_pop = pd.DataFrame(pop_data, columns=PopFitness._fields)
    with sqlite3.connect(outfile) as conn:
        df.to_sql('data', conn, if_exists='replace')
    with sqlite3.connect(outfile_pop) as conn:
        df_pop.to_sql('data', conn, if_exists='replace')
    
