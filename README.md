# insideOutside
A simple method for identifying interior and exterior points of a 3D point cloud, points. 

Method overview:
isideOutside computes a convex hull over the set of points via Delaunay triangulation. 
For each point, it then computes the set of disctances from that point to each face of the hull.
The minimum distance to the surface and variance of distances to the surface is then calculated.
Hierarchical lustering is then performed on the two dimensoinal parameter space to classify the points as either inside or outside.

input: points = mx3 matrix of cartesian coordinates. 
output: idx = mx1 bit vector indexing inside points (0) and outside points (1).
        clustData = mx2 maxtrix of parameter values. 
                    Colum 1 is the log10 minimum distance to surface. 
                    Column 2 is the variance of distances to the surface.

