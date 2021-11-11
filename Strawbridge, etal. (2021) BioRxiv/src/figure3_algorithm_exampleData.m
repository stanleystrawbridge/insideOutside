%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data examples 

set(0, 'DefaultFigureRenderer', 'painters');

src_folder = fullfile(pwd,'src');
save_folder = fullfile(pwd,'figures','figure3');

if ~exist(save_folder,'dir')
    mkdir(save_folder)
end

addpath(src_folder)
addpath(save_folder)

set(0, 'DefaultFigureRenderer', 'painters');

src_folder = fullfile(pwd,'src');
save_folder = fullfile(pwd,'figures','figure3');

if ~exist(save_folder,'dir')
    mkdir(save_folder)
end

addpath(src_folder)
addpath(save_folder)

% Set Parameters for sphere
r_sphere = 1;
O_sphere = [0 0 0];

n_sphere_points = 100;

n_spheres = 500;

% Set Parameters for ball
r_ball = [0.1 0.25 0.5 0.75 0.9 1]';
O_ball = [0 0 0];

n_ball_points = 50;

n_ball_radii = numel(r_ball);

% Set ground truth for inside outside points
ground_truth = [zeros(n_ball_points,1); %inside points
                ones(n_sphere_points,1)]; %outside points

% Make sphere
[S] = makeSphereRandom(r_sphere,n_sphere_points, O_sphere,0.0);    

axis_lims = [-1.3 1.3 -1.3 1.3 -1.3 1.3];

for i = 1:n_ball_radii
     
    % Make Ball
    [B] = makeBallRandom(r_ball(i),n_ball_points, O_ball);


    % Perform insideOutside classification
    P = [B; S];        
    [idx, clusterData] = insideOutside(P);
    
    % Set up inputs for makePlots
    Sphere = struct('r', r_sphere,'O', O_sphere);
    Ball =  struct('r', r_ball(i),'O', O_ball);
    
    obj = struct('points', P,...
                 'data',clusterData,...                 
                 'classification', idx,...                 
                 'ground_truth', ground_truth,...
                 'savePath', save_folder,...
                 'saveNameSuffix', num2str(Ball.r));
    
   makePlots(Sphere, Ball, obj, axis_lims);    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function makePlots(Sphere, Ball, obj, axis_lims)
    

    % Generate a convex hull
    DT = delaunayTriangulation(obj.points);
    vertices = convexHull(DT);

    surfaceVertices = vertices;
    surfacePoints = DT.Points;

    % Points and Shpheres    
    plotPointsAndSpheres(Sphere, Ball, obj, axis_lims)
    
    % Delaunay Triangulation and convex hull
    plotTriangulationAndHull(surfaceVertices,surfacePoints,obj, axis_lims)
    
    % Parameter Space
    plotparameterSpace(obj)
    
    % Classified points
    plotClassification(surfaceVertices,surfacePoints,obj, axis_lims)
    
end


function plotPointsAndSpheres(Sphere, Ball, obj, axis_lims)

    fig = figure(1);
    clf
    hold on

    points = obj.points;
    surface_color = [1.0 0.5 0.4];
    ball_color = [0.4 0.6 0.9];
    
    % SURFACE==============================================================
    % plot surface points--------------------------------------------------
    surface_idx = obj.ground_truth==1;
    h_outside = scatter3(points(surface_idx,1),points(surface_idx,2),points(surface_idx,3),...
        100,surface_color,'filled','MarkerFaceAlpha',0.7);
    
    % plot surface shell---------------------------------------------------
    
    [X,Y,Z] = sphere;
    
    [x_grid, y_grid] = size(Z);
    
    X = Sphere.r.*X + Sphere.O(1);
    Y = Sphere.r.*Y + Sphere.O(2);
    Z = Sphere.r.*Z + Sphere.O(3);
    
    c_idx = floor(255*(Z+ Sphere.r)/(2*Sphere.r))+1;
    n_shades = 256;
    cmap_ball =[linspace(0,surface_color(1), n_shades);
                linspace(0,surface_color(2), n_shades);
                linspace(0,surface_color(3), n_shades)]'; 
            
    CO(:,:,1) = reshape(cmap_ball(c_idx,1),x_grid, y_grid);
    CO(:,:,2) = reshape(cmap_ball(c_idx,2),x_grid, y_grid);
    CO(:,:,3) = reshape(cmap_ball(c_idx,3),x_grid, y_grid);
    
    surf(X,Y,Z, CO, 'EdgeColor','none','FaceAlpha',0.3);
    shading interp
    
    % BALL=================================================================
    % plot inner ball Points-----------------------------------------------
	ball_idx = obj.ground_truth==0;
   h_inside = scatter3(points(ball_idx,1),points(ball_idx,2),points(ball_idx,3),...
        100,ball_color,'filled','MarkerFaceAlpha',0.7);
    
    % plot inner ball shell------------------------------------------------
    if Ball.r ~= Sphere.r
        [X,Y,Z] = sphere;

        [x_grid, y_grid] = size(Z);

        X = Ball.r.*X + Ball.O(1);
        Y = Ball.r.*Y + Ball.O(2);
        Z = Ball.r.*Z + Ball.O(3);

        c_idx = floor(255*(Z+ Ball.r)/(2*Ball.r))+1;
        n_shades = 256;
        cmap_ball =[linspace(0,ball_color(1), n_shades);
                    linspace(0,ball_color(2), n_shades);
                    linspace(0,ball_color(3), n_shades)]'; 

        CO(:,:,1) = reshape(cmap_ball(c_idx,1),x_grid, y_grid);
        CO(:,:,2) = reshape(cmap_ball(c_idx,2),x_grid, y_grid);
        CO(:,:,3) = reshape(cmap_ball(c_idx,3),x_grid, y_grid);

        surf(X,Y,Z,CO,'EdgeColor','none','FaceAlpha',0.3);
        shading interp
    end
    % Make it look nice====================================================
    xlabel 'x-axis'
    ylabel 'y-axis'
    zlabel 'z-axis'

    view(3)
    grid on
    set(gca,'XTick',-1:0.5:1,'YTick', -1:0.5:1, 'ZTick', -1:0.5:1)

    set(gca,'FontSize',12,'FontName','Arial')
    
    axis equal
    axis([axis_lims])
    
    legend([h_inside,h_outside], {'Inside', 'Outside'},...
        'Orientation','horizontal','Location','north')
    
    % Saving===============================================================
    saveName = strcat('exampleData_points_and_spheres',...
                      '_ballRadius_',obj.saveNameSuffix);
    saveFigure(fig,saveName,obj.savePath) 

end


function plotTriangulationAndHull(surfaceVertices,surfacePoints,obj, axis_lims)

    fig = figure(2);
    clf
    hold on
    
    n_shades = 256;
    color = [1.0 0.5 0.4];
    cmap = [linspace(0, 1.0, n_shades);
            linspace(0, 0.5, n_shades);
            linspace(0, 0.4, n_shades)]'; 

    points = obj.points;
    
    % Plot convex hull
    trisurf(surfaceVertices,...
            surfacePoints(:,1),...
            surfacePoints(:,2),...
            surfacePoints(:,3),...
            'EdgeColor',color,'LineWidth', 3,...
            'facealpha',.3,'EdgeAlpha',0.3);

    colormap(cmap)

    % Plot all points
    scatter3(points(:,1),points(:,2),points(:,3),...
        100,'k','filled','MarkerFaceAlpha',0.7);

    xlabel 'x-axis'
    ylabel 'y-axis'
    zlabel 'z-axis'

    view(3)
    grid on
    set(gca, 'XTick', -1:0.5:1, 'YTick', -1:0.5:1, 'ZTick', -1:0.5:1)
    
    set(gca,'FontSize',12,'FontName','Arial')
    
    axis equal
    axis([axis_lims])

    
    % Saving===============================================================
    saveName = strcat('exampleData_points_and_hull',...
                      '_ballRadius_',obj.saveNameSuffix);
    saveFigure(fig,saveName,obj.savePath)    
    
end


function plotparameterSpace(obj)

    fig = figure(3);
    clf
    hold on
    
    color = [0.4 0.6 0.9;
             1.0 0.5 0.4;
             0.6 1.0 0.6;
             0.2 0.3 0.3];
    
    trueIn = (obj.classification == 0) & (obj.ground_truth == 0);
    trueOut = (obj.classification == 1) & (obj.ground_truth == 1);
    falseOut = (obj.classification == 1) & (obj.ground_truth == 0);
    falseIn = (obj.classification == 0) & (obj.ground_truth == 1);
    
    h_trueIn = scatter(obj.data(trueIn,1),obj.data(trueIn,2),100,color(1,:),...
        'filled','MarkerFaceAlpha',0.7);
    
    h_trueOut = scatter(obj.data(trueOut,1),obj.data(trueOut,2),100,color(2,:),...
        'filled','MarkerFaceAlpha',0.7);
    
    h_falseOut = scatter(obj.data(falseOut,1),obj.data(falseOut,2),100,color(3,:),...
        'filled','MarkerFaceAlpha',0.7);
    
    h_falseIn = scatter(obj.data(falseIn,1),obj.data(falseIn,2),100,color(4,:),...
        'filled','MarkerFaceAlpha',0.7);     
    
    % Make it look nice====================================================
    xlabel 'min(d)'
    ylabel 'var(d)'
    title 'Parameter Space'
    grid on

    set(gca,'XTick',0:0.25:1,'YTick', 0:0.25:1)
    
    set(gca,'FontSize',15,'FontName','Arial')
    
    axis equal
    axis([0 1 0 1])
    
    % Saving===============================================================
    saveName = strcat('exampleData_parameter_space',...
                      '_ballRadius_',obj.saveNameSuffix);
    saveFigure(fig,saveName,obj.savePath)   
    
end


function plotClassification(surfaceVertices,surfacePoints,obj, axis_lims)
    
    fig = figure(4);
    clf
    hold on
    
    n_shades = 256;
    
    color = [0.4 0.6 0.9;
             1.0 0.5 0.4;
             0.6 1.0 0.6;
             0.2 0.3 0.3];
    
    cmap = [linspace(0, color(2,1), n_shades);
            linspace(0, color(2,2), n_shades);
            linspace(0, color(2,3), n_shades)]'; 

        
    trueIn = (obj.classification == 0) & (obj.ground_truth == 0);
    trueOut = (obj.classification == 1) & (obj.ground_truth == 1);
    falseOut = (obj.classification == 1) & (obj.ground_truth == 0);
    falseIn = (obj.classification == 0) & (obj.ground_truth == 1);
        
    h_trueIn = scatter3(obj.points(trueIn,1),...
                        obj.points(trueIn,2),...
                        obj.points(trueIn,3),...
                        100,color(1,:),'filled','MarkerFaceAlpha',0.7);
                    
    h_trueOut = scatter3(obj.points(trueOut,1),...
                        obj.points(trueOut,2),...
                        obj.points(trueOut,3),...
                        100,color(2,:),'filled','MarkerFaceAlpha',0.7);
        
    h_falseOut = scatter3(obj.points(falseOut,1),...
                        obj.points(falseOut,2),...
                        obj.points(falseOut,3),...
                        100,color(3,:),'filled','MarkerFaceAlpha',0.7);
        
    h_falseIn = scatter3(obj.points(falseIn,1),...
                        obj.points(falseIn,2),...
                        obj.points(falseIn,3),...
                        100,color(4,:),'filled','MarkerFaceAlpha',0.7);

     % Plot convex hull
    trisurf(surfaceVertices,...
            surfacePoints(:,1),...
            surfacePoints(:,2),...
            surfacePoints(:,3),...
            'EdgeColor',color(2,:),'LineWidth', 3,...
            'facealpha',.2,'EdgeAlpha',0.2);

    colormap(cmap)
    
    legend([h_trueIn, h_trueOut, h_falseOut, h_falseIn],...
        {'In','Out','False Out','False In'})
    
    % Make it look nice====================================================    % Make it look nice====================================================
    xlabel 'x-axis'
    ylabel 'y-axis'
    zlabel 'z-axis'

    view(3)
    grid on
    set(gca,'XTick',-1:0.5:1,'YTick', -1:0.5:1, 'ZTick', -1:0.5:1)

    set(gca,'FontSize',12,'FontName','Arial')
    
    axis equal
    axis([axis_lims])

    
    grid on

   
    % Saving===============================================================
    saveName = strcat('exampleData_classification',...
                      '_ballRadius_',obj.saveNameSuffix);
    saveFigure(fig,saveName,obj.savePath)   
end


function saveFigure(fig,saveName,savePath)

    fullSaveName = fullfile(savePath, saveName);
    
    savefig(fig,[fullSaveName,'.fig'])
    saveas(fig,[fullSaveName,'.svg'],'svg')
    saveas(fig,[fullSaveName,'.png'],'png')

end

