To produce a simulation:
```
time ./background-sim.py -T 10 -N 10 -w 2 -L 100 -l 10 -m 0.5 -u 5e-3 -r 2.5e-3 -a .23 -b 5.34 -s bground_sim_short.selloci -A 10 -k 2 -U 0.05 -o bground_sim_short.vcf -t bground_sim_short.trees -g bground_sim_short.log
```

Now:
```py
import msprime

def stats(ts):
    print("# nodes:", ts.num_nodes)
    print("# records:", ts.num_records)
    print("# trees:", ts.num_trees)
    print("# mutations:", ts.num_mutations)
    print("sequence length:", ts.sequence_length)
    print("samples:", ts.samples())

ts = msprime.load("bground_sim_short.trees")
stats(ts)

new_ts = ts.simplify()
stats(new_ts)

```
