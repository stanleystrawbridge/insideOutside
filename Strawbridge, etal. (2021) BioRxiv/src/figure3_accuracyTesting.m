function r_ball = figure3_accuracyTesting(method,figNum)

    clc

    set(0, 'DefaultFigureRenderer', 'painters');

    save_folder = fullfile(pwd,'figures','figure3');

    if ~exist(save_folder,'dir')
        mkdir(save_folder)
    end

    addpath(save_folder)

    % Set Parameters for sphere
    r_sphere = 1;
    n_sphere_points = 100;
    O_sphere = [0 0 0];

    n_spheres = 1000;

    % Set Parameters for ball
    n_ball_radii = 100;
    r_step = 1/n_ball_radii;
    r_ball = [r_step:r_step:1]';
    n_ball_points = 50;
    O_ball = [0 0 0];

    % Set ground truth for inside outside points
    ground_truth = [zeros(n_ball_points,1); %inside points
                    ones(n_sphere_points,1)]; %outside points

    % Step over spheres and balls to generate statistics for sensitivity and
    % specificity of algorithm
    sensitivity = zeros(n_ball_radii,n_spheres);
    specificity = sensitivity;
    precision = sensitivity;

    for i = 1:n_ball_radii

        r = r_ball(i);

        disp(['processing ball of radius: ', num2str(r),'. ',...
            num2str(i), ' of ', num2str(n_ball_radii)]);

        parfor j = 1:n_spheres

            % Make Ball
            [B] = makeBallRandom(r,n_ball_points, O_ball);

            % Make sphere
            [S] = makeSphereRandom(r_sphere,n_sphere_points, O_sphere);    

            % Perform insideOutside classification
            P = [B; S];        

            switch method
                case 'insideOutside_raw'
                    idx = insideOutside_raw(P);        
                case 'insideOutside'
                    idx = insideOutside(P);      
                case 'ellipsoidalMethod'
                    idx = ellipsoidalMethod(P)
                case 'convexHullMethod'
                    idx = convexHullMethod(P);        
            end

            % Calculate Sensitivity, Specificty, and precisions
            [stats] = calculateStats(idx, ground_truth);

            sensitivity(i,j) = stats(1);
            specificity(i,j) = stats(2);
            precision(i,j) = stats(3);

        end

    end

    fig = figure(figNum);
    clf
    plotStats(r_ball,sensitivity,specificity)
    saveFigure(fig,['accuracy_', method],save_folder)

    stats = struct('sensitivity',sensitivity,...
                   'specificity',specificity,...
                   'precision',precision);
    save(fullfile(save_folder,['accuracy_', method, '_stats']), 'stats')   

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotStats(r_ball,sensitivity,specificity)

    hold on

    h_sens = plotError(r_ball,sensitivity,[1.0 0.5 0.4],'-');
    h_spec = plotError(r_ball,specificity,[0.4 0.6 0.9],':');

    legend([h_sens, h_spec],...
           {'sensitivity','specificity','precision'},...
            'Location','southwest','Orientation','vertical')
    title 'Accuracy Testing'
    xlabel 'Inner Ball Radius'
    ylabel 'Statistic'

    grid on 

    set(gca, 'XTick', 0:.25:1, 'YTick', 0:.25:1)
    set(gca, 'FontSize', 12,'FontName','Arial')
    axis equal
    axis([0 1 0 1])

end


function h = plotError(r,x,color,lineStyle)

    mu = mean(x,2);
    SEM = std(x,[],2);

    fill([r; flip(r)],[mu+SEM;flip(mu-SEM)],color,...
        'EdgeColor','none','FaceAlpha',0.5)

    h = plot(r,mu,lineStyle,'Color',color,'LineWidth',4);

end


function saveFigure(fig,saveName,savePath)

    fullSaveName = fullfile(savePath, saveName);

    savefig(fig,fullSaveName)
    saveas(fig,fullSaveName,'svg')
    saveas(fig,fullSaveName,'png')

end