# insideOutside
A simple method for classifying interior and exterior points of a 3D point cloud. 
Originally published in Stirparo and Kurowski, et al. (2021) PNAS (https://doi.org/10.1073/pnas.2008890118).
For full method, with benchmarking, please cite Strawbridge, et al. (2021) BioRxiv (https://doi.org/10.1101/2021.11.15.468285).

The method is provided as a MATLAB file insideOutside.m, using MATLAB2021a; and a Python file inside_outside.py, using Python v3.9 with NumPy v1.22.4 and SciPy v1.10.1. 

**Method overview:**
insideOutside computes a convex hull over the set of points via Delaunay triangulation. 
For each point, it then computes the set of distances from that point to each face of the hull.
The minimum distance to the surface and variance in distances to the surface is then calculated for each point.
Hierarchical clustering is then performed over the two dimensoinal parameter space to classify the points as either inside (away from the hull) or outside (near the hull).

**input**: <br /> mx3 matrix of Cartesian coordinates, where m is the number of points. 

**output**: <br />
mx1 bit vector indexing inside points (0/False) and outside points (1/True).<br />
mx2 matrix of parameter values storing information about the distance of each point to the convex hull.
