[![Build Status](https://github.com/prunett/SLAPRP.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/prunett/SLAPRP.jl/actions/workflows/CI.yml?query=branch%3Amain)

# SLAPRP.jl

Branch-Cut-and-Price algorithm for solving the Storage Location Assignment and Picker Routing Problem.

## Description

This repository implements the code from the paper "The Storage Location Assignment and Picker Routing Problem: A Generic Branch-Cut-and-Price Algorithm" authored by T. Prunet, N. Absi and D. Cattaruzza. The paper is currently under peer-review but a preprint is available on [https://arxiv.org/abs/2407.13570](https://arxiv.org/abs/2407.13570). 

## Requirements

Running with package requires `Julia >= 1.6`. The Linear programs are solved using the commercial solver `ILOG CPLEX 20.1`. Please refer to [https://github.com/jump-dev/CPLEX.jl](https://github.com/jump-dev/CPLEX.jl) to get a licence and interface CPLEX with Julia.

## Reproducing our results

To reproduce our result, you should first clone the repository. Open a Julia REPL at its root, and run the following commands:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("test/main.jl")
```

Then you can run the experiments with the command 
```julia
run(dataset, policy)
```
where dataset is the name of the dataset and policy is the routing policy. The field data set can take the values "silva" and "guo".
The field policy can take the values "exact", "return", "sshape", "midpoint", "largest". The detailed description of the datasets and routing policies is 
provided in the related publication. Note that "return" is the only policy available for the guo dataset.

For instance, run the following command to reproduce the results for the silva dataset with the return routing policy. This will create a csv file with the detailed results.

```julia
run("silva","return")
```

Beware that this command will run the algorithm on the whole testbed, which may be time consuming.

## Reference

```
@misc{prunet2024storagelocationassignmentpicker,
      title={The Storage Location Assignment and Picker Routing Problem: A Generic Branch-Cut-and-Price Algorithm}, 
      author={Thibault Prunet and Nabil Absi and Diego Cattaruzza},
      year={2024},
      eprint={2407.13570},
      archivePrefix={arXiv},
      primaryClass={cs.DM},
      url={https://arxiv.org/abs/2407.13570}, 
}
```
