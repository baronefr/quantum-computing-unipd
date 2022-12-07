# Homework 5


## Preparation

Compile with
```bash
make init
```

Execute `init.x` with the args `N (int), output file prefix (string)`.


## Exr

Compile with
```bash
make exr
```

Execute with the args `Nt (int), tmax (double), tau (double), output file prefix (string)`; for instance
```
./exr.x 200 20d0 10d0 data/psi_
```

### dataset

Since the dataset is large, you might compute it yourself with the following commands.

On $\psi_0$ (set by default `chosen_eigenstate = 1`):
```bash
./exr.x 10000 20d0 10d0 data/psi_
./exr.x 5000 20d0 10d0 data/psi_
./exr.x 4000 20d0 10d0 data/psi_
./exr.x 3000 20d0 10d0 data/psi_
./exr.x 1000 20d0 10d0 data/psi_
./exr.x 500 20d0 10d0 data/psi_
./exr.x 200 20d0 10d0 data/psi_
./exr.x 10000 20d0 1d0 data/psi_
./exr.x 10000 20d0 4d0 data/psi_
./exr.x 10000 20d0 8d0 data/psi_
./exr.x 10000 20d0 15d0 data/psi_
```

On $\psi_3$ (set `chosen_eigenstate = 3`):
```bash
./exr.x 10000 20d0 10d0 data/psi3_
```

On $\psi_{10}$ (set `chosen_eigenstate = 10`):
```bash
./exr.x 10000 20d0 10d0 data/psi3_
```