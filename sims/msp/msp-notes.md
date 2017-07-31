Will simulate:

- 12 chromosomes
- 20 x 20 grid
- 2000 smaples total (5 samples per pop)
- chromosomes all 1e8 bp

and base of

- N=10,000 per population (4 million total)
- migration rate 1e-4 per gen
- recomb rate of 2e-8 per gen per bp
- mut rate of 1e-8 per gen per bp

## Nothing much

```
SEED=$(printf "%06d" $RANDOM); mkdir nada_$SEED; time \
    ./msp-sim.py -n 12 -k 2000 -N 1e4 -w 10 -m 1e-4 -o nada_$SEED -d $SEED -u 0.0 \
        | grep -v REJECTING &> nada_$SEED/run.log &
```

## Background selection

*background_modest*:

- `N = [5000, 5000, 5000, 5000, 5000, 10000, 10000, 10000, 10000, 20000, 20000, 20000]`

```
SEED=$(printf "%06d" $RANDOM); mkdir background_modest_$SEED; time \
    ./msp-sim.py -n 12 -k 2000 -w 10 -m 1e-4 -o background_modest_$SEED -d $SEED -u 0.0 \
        -N 5000 5000 5000 5000 5000 10000 10000 10000 10000 20000 20000 20000 \
        | grep -v REJECTING &> background_modest_$SEED/run.log &
```

*background_strong*:

- `N = [1000, 2000, 2000, 5000, 5000, 5000, 10000, 10000, 10000, 10000, 20000, 20000]`

```
SEED=$(printf "%06d" $RANDOM); mkdir background_strong_$SEED; time \
    ./msp-sim.py -n 12 -k 2000 -w 10 -m 1e-4 -o background_strong_$SEED -d $SEED -u 0.0 \
        -N 1000 2000 2000 5000 5000 5000 10000 10000 10000 10000 20000 20000
        | grep -v REJECTING &> background_strong_$SEED/run.log &
```


## Local adaptation

*local_modest*:

- `migration = 1e-4 * [0.1, 0.2, 0.2, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 3.0]`

```
SEED=$(printf "%06d" $RANDOM); mkdir background_strong_$SEED; time \
    ./msp-sim.py -n 12 -k 2000 -w 10 -N 1e4 -o background_strong_$SEED -d $SEED -u 0.0 \
        -m .00001 .00002 .00002 .00005 .0001 .0001 .0001 .0001 .0001 .0002 .0002 .0003 \
        | grep -v REJECTING &> background_strong_$SEED/run.log &
```


*local_strong*:

- `migration = 1e-4 * [0.01, 0.02, 0.02, 0.1, 0.1, 0.5, 1.0, 1.0, 1.0, 4.0, 5.0, 5.0]`

## Both

*both_modest* and *both_strong*


## Recombination rate variation


## Mutation rate variation
