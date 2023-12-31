# ATOM: An Efficient Topology Adaptive Algorithm for Minor Embedding in Quantum Computing
## Introduction
Quantum annealing (QA) has emerged as a powerful technique to solve optimization problems by taking advantages of quantum physics. In QA process, a bottleneck that may prevent QA to scale up is minor embedding step in which we embed optimization problems represented by a graph, called logical graph, to Quantum Processing Unit (QPU) topology of quantum computers, represented by another graph, call hardware graph. Existing methods for minor embedding require a significant amount of running time in a large-scale graph embedding. 
To overcome this problem, in this paper, we introduce a novel notion of adaptive topology which is an expandable subgraph of the hardware graph. From that, we develop a minor embedding algorithm, namely Adaptive TOpology eMbedding (ATOM). ATOM iteratively selects a node from the logical graph, and embeds it to the adaptive topology of the hardware graph. Our experimental results show that ATOM is able to provide a feasible embedding in much smaller running time than that of the state-of-the-art without compromising the quality of resulting embedding.

## How to cite
Please cite the paper corresponding to this repository:
```
@INPROCEEDINGS{Ngo2023,
AUTHOR="Hoang Minh Ngo and Tamer Kahveci and My T. Thai",
TITLE="{ATOM:} An Efficient Topology Adaptive Algorithm for Minor Embedding in
Quantum Computing",
BOOKTITLE="2023 IEEE International Conference on Communications; Selected Areas in
Communications: Quantum Communications and Information Technology (IEEE
ICC'23 - SAC-11 QCIT Track)",
ADDRESS="Rome, Italy",
DAYS="27",
MONTH=may,
YEAR=2023
}
```
## How to run
