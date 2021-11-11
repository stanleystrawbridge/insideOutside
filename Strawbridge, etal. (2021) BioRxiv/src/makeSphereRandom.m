function [S] = makeSphereRandom(r,varargin)

    % Markus Deserno (2004) How to generate equidistributed points on the
    % surface of a sphere
    
    % Set the number of points to distribute over the surface of a sphere
    n_points = 100;
    % Set the center of the Sphere
    Ori = [0, 0, 0];
    % Set the ammount of noise in the sphere
    noise_magnitude = 0;    
    
    if numel(varargin)>0
        n_points = varargin{1};
    end
    
    if numel(varargin)>1
        Ori = varargin{2};
    end
    
    if numel(varargin)>2
         noise_magnitude = varargin{3};
    end
    
    % Set random variables
    z = r*(2*rand(n_points,1)-1);  
    
    phi = 2*pi*rand(n_points,1);
        
    x = sqrt(r.^2-z.^2).*cos(phi);
    y = sqrt(r.^2-z.^2).*sin(phi);    
    
    S = [x,y,z] + repmat(Ori,n_points,1);
    
    % Add noise%%
    s_noise = normNoise(S,n_points,noise_magnitude);
    
    S = S + s_noise;
    
end

function s_noise =  normNoise(S,n_points, noise_magnitude)

    
    noise = noise_magnitude*abs(randn(n_points,1));
    s_noise = zeros(n_points,3);

    for i = 1:n_points
        
        vec = S(i,:);
        
        normVec = vec./norm(vec);
        s_noise(i,:) = noise(i).*normVec;
        
        
    end
    
end
