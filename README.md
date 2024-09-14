# hierarchical-mpl

Hierarchical macro placement build in OpenROAD flow.

Basic flow for Hier-RTLMP 2.0  (most of the codes already in OpenROAD)

Step 1:  Leiden clustering  + rent rule -> only for standard cell clusters (see BlobPlacement)

Step 2:  regularity + connection signature based clustering -> only for macro placement (see Hier-RTLMP)

Step 3:  Mixed-size placement to determine rough location,  add boundary penalty (use DG-RePlAce)

Step 5:  Run legalization to determine the bounding box for each macro cluster

Step 6:  Determine the bus routing for the cluster level (ILP based) -> determine the pin access for each macro clustering

Step 7:  Model the pin access for each cluster as macros with fence constraints.  Then place the macro in each cluster 

## Main Procedure

