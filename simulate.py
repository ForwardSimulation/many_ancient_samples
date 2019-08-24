import fwdpy11
import numpy as np
import sys
import argparse
# import gzip
import datetime

GENOME_LENGTH = 1.0
FITNESS_SCALING = 2.0
TWONS = 1000
SIMLEN = 20


def make_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--popsize", '-N', type=int, default=None,
                          help="Population size")
    required.add_argument("--mu", '-m', type=float, default=None,
                          help="Mutation rate")
    required.add_argument('--recrate', '-r', type=float, default=None,
                          help="Recombination rate")
    required.add_argument('--seed', type=int,
                          default=None, help="RNG seed")
    required.add_argument('--outfile', '-o', type=str,
                          help="Output file name (gzipped)")
    dfe = parser.add_argument_group("DFE options")
    dfe.add_argument("--mean", "-M", type=float, default=-
                     1.0, help="Mean 2Ns of gamma DFE")
    dfe.add_argument("--shape", "-S", type=float,
                     default=1, help="Shape of gamma DFE")
    dfe.add_argument("--proportion", '-P', type=float, default=0.0,
                     help="Proportion of beneficial mutations")
    return parser


def runsim(args):
    rng = fwdpy11.GSLrng(args.seed)

    pdict = {'gvalue': fwdpy11.Multiplicative(FITNESS_SCALING),
             'rates': (0., args.mu, None),
             'nregions': [],
             'sregions': [fwdpy11.GammaS(0., GENOME_LENGTH, 1.0-args.proportion,
                                         args.mean, args.shape, FITNESS_SCALING/2.0,
                                         scaling=2*args.popsize, label=1),
                          fwdpy11.ConstantS(0, GENOME_LENGTH, args.proportion, TWONS,
                                            FITNESS_SCALING/2.0, label=2,
                                            scaling=2*args.popsize)],
             'recregions': [fwdpy11.PoissonInterval(0, GENOME_LENGTH, args.recrate)],
             'demography': np.array([args.popsize]*SIMLEN*args.popsize, dtype=np.uint32),
             # This could easily be True for these sims:
             'prune_selected': False
             }
    params = fwdpy11.ModelParams(**pdict)

    pop = fwdpy11.DiploidPopulation(args.popsize, GENOME_LENGTH)

    sampler = fwdpy11.RandomAncientSamples(
        args.seed, args.popsize, [i for i in range(10*pop.N, SIMLEN*pop.N)])

    # With a lot of ancient samples:
    # 1. RAM use already skyrockets
    # 2. Simplification slows down
    # So, we should do it a little less often:
    fwdpy11.evolvets(rng, pop, params, 1000, sampler,
                     suppress_table_indexing=True)

    return pop


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    print("Start sim:", datetime.datetime.now())
    pop = runsim(args)
    print("End sim:", datetime.datetime.now())
    print("Start binary dump:", datetime.datetime.now())
    # with gzip.open(args.outfile, 'wb') as f:
    #     pop.pickle_to_file(f)
    pop.dump_to_file(args.outfile)
    print("End binary dump:", datetime.datetime.now())
    # print("Start dump to tskit:", datetime.datetime.now())
    # ts = pop.dump_tables_to_tskit()
    # ts.dump('tskit_format.trees')
    # print("End dump to tskit:", datetime.datetime.now())
