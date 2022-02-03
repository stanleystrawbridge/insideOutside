# insideOutside
A simple method for classifying interior and exterior points of a 3D point cloud (MATLAB2021a). 
Originally published in Stirparo and Kurowski, et al. (2021) PNAS (https://doi.org/10.1073/pnas.2008890118).
For full method, with benchmarking, please cite Strawbridge, et al. (2021) BioRxiv (https://doi.org/10.1101/2021.11.15.468285).

**Method overview:**
insideOutside computes a convex hull over the set of points via Delaunay triangulation. 
For each point, it then computes the set of distances from that point to each face of the hull.
The minimum distance to the surface and variance in distances to the surface is then calculated for each point.
Hierarchical clustering is then performed over the two dimensoinal parameter space to classify the points as either inside (away from the hull) or outside (near the hull).

**input**: <br /> points = mx3 matrix of cartesian coordinates. 

**output**: <br />
**idx** = mx1 bit vector indexing inside points (0) and outside points (1).<br />
**clustData** = mx2 maxtrix of parameter values.
>> Column 1 is the log10 minimum distance to surface. <br />
>> Column 2 is the variance of distances to the surface.
