## ----eval=FALSE----------------------------------------------------------
#  library("devtools")
#  install_github("jw156605/SLICER")

## ----eval=FALSE----------------------------------------------------------
#  library(SLICER)
#  genes = select_genes(traj)
#  k = select_k(traj[,genes], kmin=5)
#  traj_lle = lle(traj[,genes], m=2, k)$Y
#  traj_graph = conn_knn_graph(traj_lle,5)
#  ends = find_extreme_cells(traj_graph, traj_lle)
#  start = 1
#  cells_ordered = cell_order(traj_graph, start)
#  branches = assign_branches(traj_graph,start)

