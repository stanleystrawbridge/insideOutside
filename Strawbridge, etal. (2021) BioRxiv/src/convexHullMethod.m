function idx = convexHullMethod(points)
    % Used by IVEN
    % Forsyth, et al. (2021) PLOS Biology
    
    idx = zeros(size(points,1),1);
    
    DT = delaunayTriangulation(points);
    vertices = convexHull(DT); 
    
    surfacePointIndices = unique(vertices(:));
    
    idx(surfacePointIndices) = 1;
    
end