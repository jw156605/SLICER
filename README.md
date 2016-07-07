SLICER 0.2.0

## Background

SLICER is an algorithm for constructing trajectories that describe gene expression changes during a sequential biological process. SLICER can capture highly nonlinear gene expression changes, automatically select genes related to the process, and detect multiple branch and loop features in the trajectory. Although the algorithm was initially developed specifically for single cell RNA-seq data, SLICER is useful for a wide range of data types including bulk RNA-seq data. 

##Installation
```{r,eval=FALSE}
library("devtools")
install_github("jw156605/SLICER")
```

## Sample Data and Code
A sample dataset containing 300 simulated "cells" each expressing 500 "genes" is included with the SLICER R package. The example below shows how to run SLICER on this sample data. Note that documentation for each function is available from within R.

```{r,eval=FALSE}
library(SLICER)
genes = select_genes(traj)
k = select_k(traj[,genes], kmin=5)
traj_lle = lle(traj[,genes], m=2, k)$Y
traj_graph = conn_knn_graph(traj_lle,5)
ends = find_extreme_cells(traj_graph, traj_lle)
start = 1
cells_ordered = cell_order(traj_graph, start)
branches = assign_branches(traj_graph,start)
```
