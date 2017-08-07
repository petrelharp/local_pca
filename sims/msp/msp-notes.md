# Timing

With migration at 1e-3, on a `10x10` grid, with `N` number per deme, recombination rate `2.5e-8` per bp, `time` in seconds, and `mem` in Kb:
```
L     N     time     mem
10e4  1000  26.10    59984
10e4  2000  92.40    60840
10e4  500   8.93     56216
10e5  1000  179.54   103048
10e5  2000  592.29   149404
10e5  500   60.58    77736
10e6  1000  5969.11  536168
2e6   1000  471.36   150432
4e6   1000  1354.48  246456
5e4   1000  13.17    58104
5e4   2000  46.03    60268
5e4   500   6.24     61340
5e5   1000  77.93    77392
5e5   2000  255.87   100880
5e5   500   27.36    67292
5e6   1000  1912.20  292120
7e4   1000  18.50    56896
7e4   2000  64.98    61832
7e4   500   7.26     53760
7e5   1000  118.36   86792
7e5   2000  365.53   118672
7e5   500   39.84    73420
8e4   1000  20.56    57240
8e4   2000  76.11    62488
8e4   500   8.24     53756
8e5   1000  135.46   91124
8e5   2000  450.49   129576
8e5   500   46.21    74836
```

This predicts that at 1000 individuals per deme (1e5 total),
doing one Morgan would require 1.89G and take 21.15 hours.



# Main simulations

Will simulate:

- 16 chromosomes
- 10 x 10 grid
- 2000 smaples total (20 samples per pop)
- chromosomes all 0.5e8 bp

and base of

- N=40,000 per population (4 million total)
- migration rate 2.5e-5 per gen
- recomb rate of 2.5e-8 per gen per bp
- mut rate of 1e-8 per gen per bp


## Nothing much

```
SEED=$(printf "%06d" $RANDOM); mkdir nada_$SEED; \
    /usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
    ./msp-sim.py -n 12 -k 2000 -N 1e4 -w 10 -m 1e-4 -o nada_$SEED -d $SEED -u 0.0 -j 16\
        | grep -v REJECTING &> nada_$SEED/run.log &
```

## Background selection

*background_modest*:

- `N = [5000, 5000, 5000, 5000, 5000, 5000, 10000, 10000, 10000, 10000, 10000, 20000, 20000, 20000, 20000, 20000]`

```
SEED=$(printf "%06d" $RANDOM); mkdir background_modest_$SEED; \
    /usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
    ./msp-sim.py -n 12 -k 2000 -w 10 -m 1e-4 -o background_modest_$SEED -d $SEED -u 0.0 -j 16\
        -N 5000 5000 5000 5000 5000 5000 10000 10000 10000 10000 10000 20000 20000 20000 20000 20000 \
        | grep -v REJECTING &> background_modest_$SEED/run.log &
```

*background_strong*:

- `N = [1000, 1000, 2000, 2000, 5000, 5000, 5000, 10000, 10000, 10000, 10000, 20000, 20000, 10000, 20000, 20000]`

```
SEED=$(printf "%06d" $RANDOM); mkdir background_strong_$SEED; \
    /usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
    ./msp-sim.py -n 12 -k 2000 -w 10 -m 1e-4 -o background_strong_$SEED -d $SEED -u 0.0 -j 16\
        -N 1000 1000 2000 2000 5000 5000 5000 10000 10000 10000 10000 10000 20000 20000 20000 20000 \
        | grep -v REJECTING &> background_strong_$SEED/run.log &
```


## Local adaptation

*local_modest*:

- `migration = 1e-4 * [0.1, 0.1, 0.2, 0.2, 0.2, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0]`

```
SEED=$(printf "%06d" $RANDOM); mkdir local_modest_$SEED; \
    /usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
    ./msp-sim.py -n 12 -k 2000 -w 10 -N 1e4 -o local_modest_$SEED -d $SEED -u 0.0 -j 16\
        -m .00001 .00001 .00002 .00002 .00002 .00005 .0001 .0001 .0001 .0001 .0001 .0002 .0002 .0002 .0003 .0003 \
        | grep -v REJECTING &> local_modest_$SEED/run.log &
```


*local_strong*:

- `migration = 1e-4 * [0.01, 0.01, 0.02, 0.02, 0.02, 0.1, 0.1, 0.5, 1.0, 1.0, 1.0, 4.0, 4.0, 5.0, 5.0, 5.0]`
```
SEED=$(printf "%06d" $RANDOM); mkdir local_strong_$SEED; \
    /usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
    ./msp-sim.py -n 12 -k 2000 -w 10 -N 1e4 -o local_strong_$SEED -d $SEED -u 0.0 -j 16\
        -m .000001 .000001 .000002 .000002 .000002 .00001 .00001 .00005 .0001 .0001 .0001 .0004 .0004 .0005 .0005 .0005 \
        | grep -v REJECTING &> local_strong_$SEED/run.log &
```


## Both

*both_modest* and *both_strong*


## Recombination rate variation


## Mutation rate variation
