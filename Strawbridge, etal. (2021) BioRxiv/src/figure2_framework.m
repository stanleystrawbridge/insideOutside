%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0, 'DefaultFigureRenderer', 'painters');


src_folder = fullfile(pwd,'src');
save_folder = fullfile(pwd,'figures','figure2');

if ~exist(save_folder,'dir')
    mkdir(save_folder)
end

addpath(src_folder)
addpath(save_folder)


% Set Parameters for surface of sphere
r_sphere = 1;
n_sphere_points = 100;
O_sphere = [0 0 0];

% Set parameters for test points
n_test_points = 50;


% Generate Surfaces and calculate distances:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Analytical solution------------------------------------------------------
test_points = [zeros(2,n_test_points);
               linspace(O_sphere(3), r_sphere + O_sphere(3),n_test_points)]';
           
d_analytical =  analyticalDistances(r_sphere, O_sphere, test_points);

% Equidistant Points on a Sphere-------------------------------------------
[S_equidistant] = makeSphereEquidistant(r_sphere,n_sphere_points, O_sphere);
test_points_equidistant = [zeros(2,n_test_points);
               linspace(O_sphere(3), r_sphere + O_sphere(3),n_test_points)]';
    
d_equidistant = calculatePhaseSpace(test_points_equidistant, S_equidistant);

% Random Uniform Points on a sphere----------------------------------------
n_samples = 1000;
D_random = zeros(n_test_points,2,n_samples);
psi = zeros(n_test_points,1);
for i = 1:n_samples
    
    [S_random] = makeSphereRandom(r_sphere,n_sphere_points, O_sphere);
    
    
    p_max = S_random(end,:);
    test_points_random = [linspace(O_sphere(1), p_max(1),n_test_points);
                   linspace(O_sphere(2), p_max(2),n_test_points);
                   linspace(O_sphere(3), p_max(3),n_test_points)]';
    
    d_random = calculatePhaseSpace(test_points_random, S_random);
    
    D_random(:,:,i) = d_random;
    psi(i) = sphericity(S_random);
    
end



% Plot Spheres%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Equidistant Points on a Sphere-------------------------------------------
fig = figure(1);
plotSurface(S_equidistant,test_points_equidistant,'Equidistant Points')
saveFigure(fig,'sphere_EquiDistant',save_folder)

% Uniform Random Points on a Sphere----------------------------------------
fig = figure(2);
plotSurface(S_random,test_points_random,'Random Points')
saveFigure(fig,'random_EquiDistant',save_folder)


% Plot Phase space%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure(3);
plotPhaseSpace(n_samples,d_analytical,D_random,d_equidistant,fig,save_folder)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d_analytical = analyticalDistances(r, Ori, test_points)
    
    points = test_points - repmat(Ori,size(test_points,1),1);
    
    x = sqrt(sum(points.^2,2));
    
    % minimudistance to surface--------------------------------------------
    min_d = r - x;

    % variance of distance to surface--------------------------------------
    
    var_d = r.^2 + x.^2 - ((3*r.^2 + x.^2)/(3.*r)).^2 ;        
    
    d_analytical = [min_d, var_d];
    
end

function [S] = makeSphereEquidistant(r,varargin)
    % Markus Deserno (2004) How to generate equidistributed points on the
    % surface of a sphere
    
    % Set the number of points to distribute over the surface of a sphere
    n_points = 1000;
    % Set the center of the Sphere
    Ori = [0, 0, 0];
    
    if numel(varargin)>0
        n_points = varargin{1};
    end
    
    if numel(varargin)>1
        Ori = varargin{2};
    end
    
    S = zeros(n_points,3);
    
    N_count = 0;

    a = (4*pi*(r.^2))/n_points;
    d = sqrt(a);

    M_theta = ceil(pi./d);

    d_theta = pi./M_theta;
    d_phi = a./d_theta;
            
    for m = 1:M_theta
        
        theta = pi*(m+0.5)/M_theta;
        M_phi = ceil(2*pi*sin(theta)./d_phi);
        
        for n = 1:M_phi
        
            phi = 2*pi*n/M_phi;
            
            N_count = N_count + 1;
            
            S(N_count,:) = r.*[sin(theta).*cos(phi),...
                               sin(theta).*sin(phi),...
                               cos(theta)];
                        
        end
        
    end  

    S = S(1:N_count,:) + repmat(Ori,N_count,1);   
    
    p_cap = [0 0 r]+ [0 0 Ori(3)];
    
    if ~any(all(ismember(S,p_cap),2))
        S = [S;p_cap];
    end
    
end


function [S] = makeSphereRandom(r,varargin)

    % Markus Deserno (2004) How to generate equidistributed points on the
    % surface of a sphere
    
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
    
    % Set random variables
    z = r*(2*rand(n_points,1)-1);  
    
    phi = 2*pi*rand(n_points,1);
    
    
    x = sqrt(r.^2-z.^2).*cos(phi);
    y = sqrt(r.^2-z.^2).*sin(phi);
    
    S = [x,y,z] + repmat(Ori,n_points,1);
    
end


function [d] = calculatePhaseSpace(P,S) 

    n_points = size(P,1); 
    
    d = zeros(n_points, 2);
    
    for i = 1:n_points
       
        [d_min, var_d] = distance2Surf(P(i,:),S);
        
        d(i,:) = [d_min, var_d];
        
    end

end


function [d_min, var_d] = distance2Surf(p,S)
    
    n_points = size(S,1);
    
    d = sqrt(sum((S - repmat(p,n_points,1)).^2,2));
    
    d_min = min(d);
    var_d = var(d);
    
end

function psi = sphericity(S)

    shp = alphaShape(S,inf);
    
    V = shp.volume;
    A = shp.surfaceArea;
    
    psi = ((pi^(1/3)).*((6*V).^(2/3)))./A;

end


function plotSurface(S,test_points,T)

    clf

    alpha = 0.5;

    subplot(121)
    hold on
    scatter3(S(:,1),S(:,2),S(:,3),75,...
        'filled','MarkerFaceAlpha',alpha)
    
    plotTestPoint(test_points);

    [X,Y,Z] = sphere;
    surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.3)
    shading interp
    
    makeItLookNice(T,3)
    
    legend('Surface Points','Test Points','Origin',...
        'Location','southoutside')
    
    set(gca,'FontName','Arial')    
    
    % Plot Orthogonal views
    angle = [  0.0, 90;
               0.0,  0;
              90.0 , 0];
          
    for i = 1:3
        
        subplot(3,2,2*i)
        hold on
        scatter3(S(:,1),S(:,2),S(:,3),75,...
            'filled','MarkerFaceAlpha',alpha)
        
        plotTestPoint(test_points);
            
        [X,Y,Z] = sphere;
        surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.3)
        shading interp
        
        makeItLookNice('', angle(i,:))
        
        set(gca,'FontName','Arial')    
        
    end
    
    subplot(3,2,2)
    title('Orthogonal Views')
   
    set(gcf,'renderer','painters');

end


function plotTestPoint(test_points)

    scatter3(test_points(:,1),test_points(:,2),test_points(:,3),25,'r',...
        'filled','MarkerFaceAlpha',0.5)
     
    scatter3(test_points(1,1),test_points(1,2),test_points(1,3),50,'k','filled')
    
end


function makeItLookNice(T,v)

    axis equal
    view(v)
    xlabel x-axis
    ylabel y-axis
    zlabel z-axis
    title(T)
    grid on
    
    axis([-1 1 -1 1 -1 1])
    
    set(gca,'XTick',[-1 0 1])
    set(gca,'YTick',[-1 0 1])
    set(gca,'ZTick',[-1 0 1])
    
end


function saveFigure(fig,saveName,savePath)

    fullSaveName = fullfile(savePath, saveName);
    
    savefig(fig,fullSaveName)
    saveas(fig,fullSaveName,'svg')
    saveas(fig,fullSaveName,'png')

end


function plotPhaseSpace(n_samples,d_analytical,D_random,d_equidistant,fig,save_folder)

    clf
    hold on

    for i = 1:n_samples 
        plot(D_random(:,1,i),D_random(:,2,i),'color',[0.5, 0.5, 0.5, 10/n_samples],'LineWidth',4)
    end

    h_analytical = plot(d_analytical(:,1),d_analytical(:,2),'r','LineWidth',4);
    h_equi = plot(d_equidistant(:,1),d_equidistant(:,2),'k:','LineWidth',4);
    
    h_rand = plot(0,0,'color',[0.5, 0.5, 0.5],'LineWidth',4);
    shg

    plot([0 1], (2/9)*[1 1],'k--','LineWidth',1.5);
    xlabel('Minimum Distance to Surface')
    ylabel('Variance of Distances to Surface')
    title('Inverse Relationship in Phase Space')

    legend([h_analytical,h_equi,h_rand],{'Analytical','Equidistant','Random'})
    set(gca,'FontSize',15,'FontName','Arial')
    set(gca,'XTick', 0:0.25:1)
    set(gca,'YTick', 0:0.1:.4)
  
    axis([0 1 0 .4])
    
    set(gcf,'renderer','painters');
    
    saveFigure(fig,'phaseSpace',save_folder)

end
