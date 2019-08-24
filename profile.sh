#!/bin/bash

# This script depends on paths on KRT's laptop,
# or any Ubuntu/Pop OS 19.04 installation, presumably.

# Requires fwdpy11 to be built with profiling enabled:
# python3 setup.py build_ext -i --enable-profiling
# Requires google perftools to be installed:
# sudo apt install gperftools

# Clean out any existing profiles:
rm -f profile*.pdf

PYTHONPATH=$HOME/src/fwdpy11 LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libprofiler.so CPUPROFILE=simulate.prof python3 \
    simulate.py --popsize 1000 --mu 1e-3 --recrate 0.5 --mean -1 --proportion 1e-3 --seed 101 --outfile pop.bin

~/go/bin/pprof --pdf `which python3` simulate.prof

PYTHONPATH=$HOME/src/fwdpy11 LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libprofiler.so CPUPROFILE=freqchange.prof python3 \
    freq_change.py pop.bin freqs.sqlite3 0.0 1254636

~/go/bin/pprof --pdf `which python3` freqchange.prof
