SLICER 0.2.0

INSTALLATION:

library("devtools")

install_github("jw156605/SLICER")

USAGE:

Documentation for the SLICER package, including example commands for analyzing the sample data distributed with SLICER, is available from within R. To access the documentation for a function within the package, type help(function_name) or ?function_name. For example, ?select_k will load the documentation for the select_k function. The sample data bundled with SLICER is a simulated trajectory that includes a branch and is stored in a matrix named "traj". The example code in the documentation for each function refers to this dataset.

The code segment below shows how to run SLICER on the sample data.

library(SLICER)
genes = select_genes(traj)
k = select_k(traj[,genes], kmin=5)
traj_lle = lle(traj[,genes], m=2, k)$Y
traj_graph = conn_knn_graph(traj_lle,5)
ends = find_extreme_cells(traj_graph, traj_lle)
start = 1
cells_ordered = cell_order(traj_graph, start)
branches = assign_branches(traj_graph,start)

