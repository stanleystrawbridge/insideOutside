clc

set(0, 'DefaultFigureRenderer', 'painters');

save_folder = fullfile(pwd,'figures','figure5');

if ~exist(save_folder,'dir')
    mkdir(save_folder)
end

addpath(save_folder)

data_files = dir(fullfile(save_folder,'*stats.mat'));
n_data_files = numel(data_files);

% plot and save classification rates

colors = [1.0 0.5 0.4 ; % Salmon#
          0.4 0.6 0.9]; % Cornflower Blue       
      
      
salmonMap = [linspace(0,1,256);
             linspace(0,.5,256);
             linspace(0,.4,256)]';

cornflowerMap = [linspace(0,.4,256);
                 linspace(0,.6,256);
                 linspace(0,.9,256)]'; 
%%
for i = 1:n_data_files
    
    method = extractBetween(data_files(i).name,'accuracy_','_stats.mat');
    save_path = fullfile(save_folder,method{:});
    
    load(fullfile(data_files(i).folder,data_files(i).name));
    
    save_name = [method{:}, '_outsideAccuracy'];    
    plotMap(stats.r_ball,stats.noise,stats.sensitivity,salmonMap,save_path,save_name,[0 0.3]);   
    
    save_name = [method{:}, '_insideAccuracy'];
    plotMap(stats.r_ball,stats.noise,stats.specificity,cornflowerMap,save_path,save_name,[0 0.15]);
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plotMap(x,y,z,map,save_folder,save_name,sigmaClims)
    
    % Plot mean------------------------------------------------------------
    fig = figure(1);
    clf
    hold on
    
    mu = mean(z,3)';
    h = pcolor(x,y,mu);   
    set(h, 'EdgeColor', 'none');
    caxis([0,1])
    colormap(map)
    colorbar    

    [C,h] = contour(x,y,mu,[0:0.1:0.9],'LineColor','k', 'LineWidth',2);
    clabel(C,h,'FontSize',20);
        
    
    axis square
    
    xlabel('Inner Ball Radius')
    ylabel('Noise Factor')    

    set(gca,'FontName','Arial','FontSize',12)
    
    saveFigure(fig,save_folder,[save_name, '_mu'])
    
    % Plot Std-------------------------------------------------------------
    fig = figure(2);
    clf    
    hold on
    
    sigma = std(z,[],3)';
    h = pcolor(x,y,sigma);
    set(h, 'EdgeColor', 'none');   
    clims = caxis;
    caxis(sigmaClims)
    colormap(map)
    colorbar
        
%     contour(x,y,sigma,'LineColor','k', 'LineWidth',2)
    
    axis square
    
    xlabel('Inner Ball Radius')
    ylabel('Noise Factor')    

    set(gca,'FontName','Arial','FontSize',12)
    
    saveFigure(fig,save_folder,[save_name, '_sigma'])
    
end

function saveFigure(fig,savePath,saveName)

    if ~exist(savePath,'dir')
        mkdir(savePath)
    end

    fullSaveName = fullfile(savePath, saveName);
    
    savefig(fig,fullSaveName)
    saveas(fig,fullSaveName,'svg')
    saveas(fig,fullSaveName,'png')

end
