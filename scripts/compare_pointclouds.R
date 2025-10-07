#--------------------------------------------------------------------------
# Name:         compare_pointcloud.R
# Description:  comparing the harmonized laz files in the point cloud level
#               checking pulse count etc. 
# Author:       Svenja Dobelmann
#--------------------------------------------------------------------------

# TBC! 
# compare point clouds

library(lasR)
library(lidR)


folder = "R:/AG_Magdon/datensaetze/solling/dobelmann/leaf-on_leaf-off_data/04_harmonized/lon23"
ctg = readLAScatalog(folder)


pipeline <- summarise(metrics = "count")


ans = exec(pipeline, on = ctg, progress = T)

