function [B] = makeBallRandom(r,varargin)

    
    % Set the number of points to distribute over the surface of a sphere
    n_points = 100;
    % Set the center of the Sphere
    Ori = [0, 0, 0];
    
    if numel(varargin)>0
        n_points = varargin{1};
    end
    
    if numel(varargin)>1
        Ori = varargin{2};
    end
    
    B = zeros(n_points,3); 
    count = 1;
    while count <= n_points
        
        p = r*(2*rand(1,3)-1);
        r_p = sqrt(sum(p.^2));
    
        if r_p < r
            
           B(count,:) = p;
            count = count + 1;
            
        end
        
    end
   
    B = B - repmat(Ori,n_points,1);
    
end