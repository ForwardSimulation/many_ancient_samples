import fwdpy11
import sys
import gzip
import pandas as pd
import sqlite3
from collections import namedtuple

AlleleFreq = namedtuple(
    'AlleleFreq', ['generation', 'position', 'origin', 'dac',
                   'fixation', 'neutral', 'label'])


def count_frequencies(pop, num_neutral):
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
    # Iterate over all samaple nodes.
    # False means to exclude the "alive" individuals
    # at time pop.generation
    for t, n, m in pop.sample_timepoints(False):
        tables, idmap = fwdpy11.simplify_tables(pop.tables, n)
        muts_visited = 0
        # NOTE: do not update samples lists below nodes.
        for tree in fwdpy11.TreeIterator(tables, idmap[n]):
            for m in tree.mutations():
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

                # Edge case: we skip mutations that fixed prior
                # to any of these time points
                fixed_before = False
                if dac == 2*pop.N and fixed == 0:
                    fixed_before = True
                if not fixed_before:
                    data.append(AlleleFreq(t, pos, origin,
                                           dac, fixed, int(m.neutral),
                                           pop.mutations[m.key].label))
        assert muts_visited == len(tables.mutations)
    return data


if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]
    neutral_mu = float(sys.argv[3])
    seed = int(sys.argv[4])

    pop = fwdpy11.DiploidPopulation.load_from_file(infile)

    rng = fwdpy11.GSLrng(seed)
    print(len(pop.mcounts))
    # NOTE: this function initiates mutation counting
    # via tree sequences as a side-effect. :(
    n = fwdpy11.infinite_sites(rng, pop, neutral_mu)
    print(len(pop.mcounts))
    print(f"{n} neutral variants added")

    freqs = count_frequencies(pop, n)
    df = pd.DataFrame(freqs, columns=AlleleFreq._fields)
    with sqlite3.connect(outfile) as conn:
        df.to_sql('data', conn, if_exists='replace')
