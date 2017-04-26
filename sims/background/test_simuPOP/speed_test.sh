#!/bin/bash

# note %M is:
#   M      Maximum resident set size of the process during its lifetime, in Kilobytes.


/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
    python3 simple-background-sim.py --generations 1000 --popsize 100 --nsamples 10 \
        --nloci 2 --recomb_rate 1e-6 --length 1e6 --ancestor_age 10 --mut_rate 0 --sel_mut_rate 0 -s /dev/stdout -e /dev/stdout &> speed_test.log


/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
    python3 simple-background-sim_no-msprime.py --generations 1000 --popsize 100 --nsamples 10 \
        --nloci 2 --recomb_rate 1e-6 --length 1e6 --ancestor_age 10 --mut_rate 0 --sel_mut_rate 0 -s /dev/stdout -e /dev/stdout &> speed_test_no_msprime.log


/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
    python3 simple-background-sim_test_output.py --generations 1000 --popsize 100 --nsamples 10 \
        --nloci 2 --recomb_rate 1e-6 --length 1e6 --ancestor_age 10 --mut_rate 0 --sel_mut_rate 0 -s /dev/stdout -e /dev/stdout -t test.recomb &> speed_test_file_output.log
