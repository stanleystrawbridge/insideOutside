function [idx, clustData] = insideOutside(points)

    DT = delaunayTriangulation(points);
    vertices = convexHull(DT);        
        
    % Calculate origin
    Ori = mean(points);
    
    % find maximum distance from Ori
    r_empirical = max(sqrt(sum((points- repmat(Ori,size(points,1),1)).^2,2)));
    
    % Calculate the distance of each cell from each face of the
    % triangulated embryo surface.   
    [~,minDistances] = distance2surf(vertices, points);

    minDist2Surf = min(abs(minDistances),[],3);
    varDist2Suf = var(abs(minDistances),[],3);

    % Cluster to get inside/outside index
    clustData = [log10(minDist2Surf+0.01);...
                 varDist2Suf]';
     
    clustData = clustData-min(clustData);
    clustData = clustData./repmat(max(clustData),size(points,1),1);
    
    % Perform Clustering
    clustIdx = clusterdata(clustData,'maxclust',2,'linkage','ward');
    cntrd = grpstats(clustData(:,1), clustIdx,'mean');
    outsideIdx = find(cntrd == min(cntrd));
    
    idx = clustIdx == outsideIdx;

end


function [projX0,minDistances] = distance2surf(vertices, points)

    % Get the number of faces and the number cells
    nPoints = size(points,1);
    nFaces = size(vertices,1);
    
    % Build array of faces
    faces = cell(1,nFaces);    
    for i = 1:nFaces
        faces{i} = points(vertices(i,:),:);
    end
    
    % Calcualte the distace of each cell to each face and store the data in
    % distances -----------------------------------------------------------
    %
    % The distance of a point (x0) from a finite plane S, given three
    % points (x1,x2,x3) in the plane is calculated by:   
    % (1) Calculate the normal vector to the plane,
    %    n = cross(x2-x1,x3-x1)/norm(cross(x2-x1,x3-x1)).
    % (2) Project the difference of x0 and one of the point onto the normal
    % vector,
    %    distance = dot(n,x0-x1).
    % (3) Determine if the proj(x0) onto the coplanar surface of S is
    % within the domain of S
    % (4) If not, Find the closet point of the bounary of S.

    % Calculate the normal vector to each faces
    normVecs = cellfun(@(f) cross(f(2,:)-f(1,:),f(3,:)-f(1,:))/...
        norm( cross(f(2,:)-f(1,:),f(3,:)-f(1,:))),faces,...
        'UniformOutput',false);
    
    % Project the difference vector onto the normal vector
    X0 = repmat(points',1,1,nFaces);
    X1 = repmat(points(vertices(:,1),:)',1,1,nPoints);
    X1 = permute(X1, [1 3 2]);
    diffMat = X0-X1;

    normMat =  repmat(cat(1,normVecs{:})',1,1,nPoints);
    normMat = permute(normMat, [1 3 2]);
    minSurfDistances = dot(normMat,diffMat);
    
    % Determine if the proj(X0) is in the triangluar domain----------------
    
    % calculate proj(X0) onto finite plane S
    surfProjX0 = X0 - normMat.*repmat(minSurfDistances,3,1,1);
    
    X2 = repmat(points(vertices(:,2),:)',1,1,nPoints);
    X2 = permute(X2, [1 3 2]);

    X3 = repmat(points(vertices(:,3),:)',1,1,nPoints);
    X3 = permute(X3, [1 3 2]);

    u = X2 - X1;
    v = X3 - X1;
    w = X0 - X1;
    
    n = cross(u,v);
    
    gamma = dot(cross(u,w),n)./dot(n,n);
    beta  = dot(cross(w,v),n)./dot(n,n);
    alpha = 1 - gamma - beta;
    
    T =   permute(gamma,[ 2 3 1])<=1 & permute(gamma,[ 2 3 1])>=0 ...
        & permute(beta, [ 2 3 1])<=1 & permute(beta, [ 2 3 1])>=0 ...
        & permute(alpha,[ 2 3 1])<=1 & permute(alpha,[ 2 3 1])>=0;
    
    surfProjX0(:,~T) = 0;
    minSurfDistances(:,~T) = 0;
    
    % Calculate the projection on the point to the line--------------------
    [p12, d12] = minDistPointOnSegment(X1,X2,X0);
    [p13, d13] = minDistPointOnSegment(X1,X3,X0);     
    [p23, d23] = minDistPointOnSegment(X2,X3,X0);
      
        
    D(:,:,1) = permute(d12,[2 3 1]);
    D(:,:,2) = permute(d13,[2 3 1]);
    D(:,:,3) = permute(d23,[2 3 1]);
    
    [minEdgeDistance,minEdgeIndex] = min(D,[],3);
    
    minEdgeDistance = permute(minEdgeDistance, [3 1 2] );
    
    minEdgeIndex = permute(minEdgeIndex, [3 1 2] );  
    
    p12(:,minEdgeIndex~=1) = 0;
    p13(:,minEdgeIndex~=2) = 0;
    p23(:,minEdgeIndex~=3) = 0;
    
    edgeProjX0 = p12+p13+p23;
    
    edgeProjX0(:,T) = 0;
    minEdgeDistance(:,T) = 0;
    
    projX0 = surfProjX0 + edgeProjX0;
    minDistances = minSurfDistances + minEdgeDistance;

end


function [point, distance] = minDistPointOnSegment(a,b,p)
    
    S = sqrt(sum((a-b).^2,1));
    A = sqrt(sum((a-p).^2,1));
    B = sqrt(sum((b-p).^2,1));
    
    % Calculate projection distance and vector
    ba = b - a; 
    pa = p - a;
    pb = p - b;
    
    dProjP = sqrt(sum(cross(ba,pa).^2,1))./ sqrt(sum(ba.^2,1));
    projP = a + repmat(dot(pa,ba)./dot(ba,ba),3,1,1).*ba;
       
    % Determing if the projected point is within the line segment
    onOrOff = A.^2 + B.^2 - 2.*(dProjP.^2) < S.^2;
    
    projP(:,~onOrOff) = 0;
    dProjP(:,~onOrOff) = 0;
    % get the minEndPoint and Distance
    endPointDistance(:,:,1) = permute(sqrt(sum(pa.^2,1)),[2 3 1]);
    endPointDistance(:,:,2) = permute(sqrt(sum(pb.^2,1)),[2 3 1]);
    
    [minEndPointDistance,minEndPointIndex] =  min(endPointDistance,[],3);
    
    minEndPointDistance = permute(minEndPointDistance, [3 1 2] );
    minEndPointIndex = permute(minEndPointIndex, [3 1 2] );
    
    a(repmat(minEndPointIndex==2,3,1,1))=0; 
    b(repmat(minEndPointIndex==1,3,1,1))=0;
    
    minEndPoint = a+b;
    minEndPoint(:,onOrOff) = 0;
    
    minEndPointDistance(:,onOrOff) = 0;
    
    % Now us onOrOff to get the correct minimem distance to the line
    % segment
    point = projP +  minEndPoint;
    distance = dProjP+ minEndPointDistance;

end