# Pinocchio_heffte



## Compilation

* [GPU compilation (enabled by default)]

```
make
```

* [CPU compilation (Pure MPI)]

```
make OMP=NO FULL_GPU=NO
```

* [CPU compilation (Hybrid MPI+OpenMP)]

```
make OMP=YES FULL_GPU=NO
```

* [CPU compilation with gcc/clang]

In the Makefile, replace this line
``` 
OMP_FLAG = -mp=multicore
```
With this one 
```
OMP_FLAG = -fopenmp
```

If your machine is neither LeonardoBoost nor magellanus, please modify SYSTYPE and its corresponding dependencies in the Makefile


