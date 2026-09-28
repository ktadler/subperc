# Spaced deletions and bond percolation

Code by Katerina Adler for the computations in *Dismantling bond-percolation subgraphs of random regular graphs*.

## The simulator

`subperc_sim` builds each instance and runs the scan in the simulations section. One instance is a Hamilton cycle on `n` vertices together with a uniform random `(d-2)`-regular graph on the same vertices. The shortcut graph is a configuration model. Self-loops, repeated pairs, and edges that repeat the cycle are removed by random double-edge swaps.

Deletions follow the floating-range rule: a counter steps by `1/D`, and the vertex `floor(counter)` is removed. The remaining edges are occupied independently with probability `p`.

```bash
c++ -O2 -std=c++17 -o subperc_sim subperc_sim.cpp
./subperc_sim 1000000 3 3 12345 out
./subperc_sim 1000000 3 3 12345 out --gap 0.1
```

Arguments are `n d nruns seed outdir`. Run `r` (counting from 0) uses the generator `seed + 1000003*r + d`. The batches in the paper are `n = 10^6`, `nruns = 3`, `seed = 12345`, and `d` in `{3, 4, 5, 6, 10}`.

Without `--gap`, the program writes `table<d>reg.tab`. The columns `p13_*` and `p23_*` are the bracket midpoints where the largest occupied component crosses `n^{1/3}` and `n^{2/3}`. At `n = 10^6` the second cutoff is `0.01n`, which is the threshold in the paper. The markers are the mean of the `p23_*` columns.

With `--gap 0.1` the program does not scan those cutoffs. For each `D` it solves `E(D,p) = 0.9` and `E(D,p) = 1.1` and writes the largest occupied component to `size<d>_gap0.1.tab`.

## The figures

`plot_compare.py` reads one directory that contains both kinds of table and writes `pvsD.png` and `size_gap_random.png`.

```bash
python3 -m pip install -r requirements.txt
python3 plot_compare.py data
```

The directory `data/` holds the tables from which those two figures were drawn.
