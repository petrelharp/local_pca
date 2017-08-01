# Timing

```

X="_parallel"
L=1e4
W=10
for K in 100 400 1000
do
  (S=$RANDOM; O=test${X}_$S; mkdir $O;
   (echo "L=$L; K=$K; W=$W; ./msp-sim.py -n 16 -k $K -N 1e4 -w $W -m 1e-4 -o $O -d $S -u 0.0 -L $L";
    /usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
       ./msp-sim.py -n 12 -k $K -N 1e4 -w 10 -m 1e-4 -o $O -d $S -u 0.0 -L $L -j 4) \
           &>> $O/timing_runs${X}.log ) &
done

(for x in test_*/timing_runs${X}.log; do 
    cat $x | tr '\n' ' ' | awk '{ print $1,$2,$3,$24}' | sed -e 's/;//g'; done) | sort -k 4 >> timing_runs$X.log

```


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
