PYTHONPATH=$HOME/src/fwdpy11 python3 simulate.py -N 250 -m 0.001 -r 0.001 --seed 42 -o test.ts
PYTHONPATH=$HOME/src/fwdpy11 python3 freq_change.py test.ts test.sql test_pop.sql 0.001 50


