%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inside cells are indexed as 0 for all methods: Sox2 GMM, MINS, IVEN, and
% insideOutside.

set(0, 'DefaultFigureRenderer', 'painters');

src_folder = fullfile(pwd,'src');
results_folder = fullfile(pwd,'figures','figure4','results');
if ~exist(results_folder,'dir')
    mkdir(results_folder)
end

data_file = fullfile(pwd,'figures','figure4','data','Oct4_deletable_sox2_staining.xlsx');
data = readtable(data_file);

colors = [1.0 0.5 0.4 ; % Salmon#
          0.4 0.6 0.9 ; % Cornflower Blue                          
          0.6 1.0 0.6 ; % Mint Green
          0.8 0.3 0.0 ; % Burnt Orange
          0.2 0.3 0.3]; %Charcoal

%% Set up data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subset data to remove overgrown embryo
data = data(data.in_out==1,:);

% Build metadata
[values, embryoSize] = grpstats([data.embryo_id, data.stage, data.genotype],data.embryo_id,{'mean','numel'});

metaData = struct('stages',unique(data.stage),...
             'genotypes',unique(data.genotype),...
             'markerNames',{{'DAPI','Oct4','Sox2'}},...
             'embryoID',values(:,1),...
             'embryoStage',values(:,2),...
             'embryoGenotype',values(:,3),...
             'embryoSize',embryoSize(:,1),...
             'resultsFolder',results_folder);

% Get fluorescence values
fluorescence = log10([data.dapi_sum, data.oct4_sum, data.sox2_sum]./repmat(data.volume,1,3)+1);
fluorescence = fluorescence-repmat(min(fluorescence),size(fluorescence,1),1);
fluorescence = fluorescence./repmat(max(fluorescence),size(fluorescence,1),1);
%%
% Set threshold using 2-component gaussian mixture modelling
fig1 = figure(100);
clf
numComponents = 2;

Sox2_index = performGaussianMixtureModelling(fluorescence(:,[1 end]),numComponents);
plotGMM(fluorescence,Sox2_index,colors,metaData)


%% Perform insideOutside Sorting algorithms%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set index vextors
% Data output from MINS
MINS_output_index = data.te_icm==2;

% initalizing other indexing vectors
ellipsoidal_index = zeros(size(data.cell_id)); 
convexHull_index = zeros(size(ellipsoidal_index ));
insideOutside_index = zeros(size(ellipsoidal_index ));
insideOutside_clusterData =zeros([size(ellipsoidal_index,1),2]);

% Step over embryos.
for i = 1:numel(metaData.embryoID)

    embryoID = data.embryo_id == metaData.embryoID(i);
    
    x = data.x(embryoID);
    y = data.y(embryoID);
    z = data.z(embryoID);       
    
    % Sox2 GMM ------------------------------------------------------------    
    fig_sox2 = figure(10);
    plotEmbryo(x,y,z,Sox2_index(embryoID),Sox2_index(embryoID),i,metaData,'sox2_gmm',fig_sox2)
        
    % MINS output ---------------------------------------------------------    
    fig_MINS = figure(11);
    plotEmbryo(x,y,z,MINS_output_index(embryoID),Sox2_index(embryoID),i,metaData,'mins',fig_MINS)
        
    % ellipsoidal method --------------------------------------------------    
    ellipsoidal_index(embryoID) = ellipsoidalMethod([x,y,z]);     
    
    fig_ellipsoidal = figure(12);
    plotEmbryo(x,y,z,ellipsoidal_index(embryoID),Sox2_index(embryoID),i,metaData,'ellipsoidal',fig_ellipsoidal)
    
    % convex hull methodf method ------------------------------------------    
    convexHull_index(embryoID) = convexHullMethod([x,y,z]);
    
    fig_convexHull = figure(13);    
    plotEmbryo(x,y,z,convexHull_index(embryoID),Sox2_index(embryoID),i,metaData,'convexHull',fig_convexHull)
    
    % insideOutside -------------------------------------------------------
    [insideOutside_index(embryoID), insideOutside_clusterData(embryoID,:)] = insideOutside([x,y,z]); 
    
    fig_insideOutside = figure(14);
    plotEmbryo(x,y,z,insideOutside_index(embryoID),Sox2_index(embryoID),i,metaData,'insideOutside',fig_insideOutside)
    
    fig_insideOutside_clustData = figure(15);
    plotInsideOutsideClusterData(insideOutside_clusterData(embryoID,:),insideOutside_index(embryoID),Sox2_index(embryoID),i,metaData,'insideOutside_clusterData',fig_insideOutside_clustData)
    
end

%%

results = struct('ground_truth',Sox2_index,...
                 'classificationMethods',{{'mins','ellipsoidal','convexHull','insideOutside'}},...
                 'classification',[MINS_output_index,ellipsoidal_index,convexHull_index,insideOutside_index],...
                 'embryo',data.embryo_id,...
                 'stage', data.stage,...
                 'genotype',data.genotype,...
                 'saveaPath',fullfile(metaData.resultsFolder,'stats'));
 
if ~exist(results.saveaPath,'dir')
    mkdir(results.saveaPath)
end

%%
 doAnalysis(results)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function groupsByGMM = performGaussianMixtureModelling(input,numComponents)

    % set GMM parameters            
    number_of_fitting_replicates = 20;
    max_itterations = 1000;
    opt = statset('MaxIter',max_itterations);

    % Sort into Sox2 high and low;
    gmm = fitgmdist(input,numComponents,...
                    'Replicates',number_of_fitting_replicates ,...
                    'Options',opt);
    
    clust = gmm.cluster(input);
    [m,idx] = sort(gmm.mu(:,2));
        
    groupsByGMM = (clust == idx(1));
                 
end


function plotGMM(fluorescence,Sox2_index,colors,metaData)

    % build color matrix
    c = zeros(size(Sox2_index,1),3);
    c(Sox2_index==0,:) = repmat(colors(2,:),sum(Sox2_index==0),1);
    c(Sox2_index==1,:) = repmat(colors(1,:),sum(Sox2_index==1),1);
    
    scatter(fluorescence(:,3),fluorescence(:,1),[],c,...
        'filled','MarkerFaceAlpha',0.5)
       
    title('GMM of Sox2 Signal')
    
    xlabel(metaData.markerNames{3})
    ylabel(metaData.markerNames{1})
 
 
    grid on
    
    set(gca,'FontSize',12,...
        'XTick',[0:0.25:1],...
        'YTick',[0:0.25:1],...
        'ZTick',[0:0.25:1])

end

function plotEmbryo(x,y,z,index,ground_truth,embryoNumber,metaData,saveName,fig)

    X = [x,y,z];
    X = X - repmat(min(X),size(X,1),1);
    
    x = X(:,1);
    y = X(:,2);
    z = X(:,3);
    colors = [1.0 0.5 0.4 ; % Salmon
              0.4 0.6 0.9 ; % Cornflower Blue                          
              0.6 1.0 0.6 ; % Mint Green
              0.2 0.3 0.3]; %Charcoal
          

    % build color matrix
    c = zeros(size(index,1),3);
    true_in = (index == 0) & (ground_truth == 0); 
    false_in = (index == 0) & (ground_truth == 1); 
    true_out = (index == 1) & (ground_truth == 1); 
    false_out = (index == 1) & (ground_truth == 0); 
    
    
    c(true_in,:) = repmat(colors(2,:),sum(true_in==1),1);
    c(false_in,:) = repmat(colors(4,:),sum(false_in==1),1);
    c(true_out,:) = repmat(colors(1,:),sum(true_out==1),1);
    c(false_out,:) = repmat(colors(3,:),sum(false_out==1),1);
        
    scatter3(x,y,z,150,c,'filled','MarkerFaceAlpha',0.7)
       
    title(['Embryo ', num2str(metaData.embryoID(embryoNumber))])
    
    xlabel x-axis
    ylabel y-axis
    zlabel z-axis
 
    grid on
    
    set(gca,'FontSize',12)
    
    embryoName = ['embryo_', sprintf('%03d', metaData.embryoID(embryoNumber))];
    
    savePath = fullfile(metaData.resultsFolder, embryoName);

    if ~exist(savePath,'dir')
        mkdir(savePath)
    end

    saveFigure(fig,saveName,savePath)
        
end

function plotInsideOutsideClusterData(clusterData,index ,ground_truth,embryoNumber,metaData,saveName,fig)
          
    colors = [1.0 0.5 0.4 ; % Salmon#
              0.4 0.6 0.9 ; % Cornflower Blue                          
              0.6 1.0 0.6 ; % Mint Green
              0.2 0.3 0.3]; %Charcoal
          
    % build color matrix
    c = zeros(size(index,1),3);
    true_in = (index == 0) & (ground_truth == 0); 
    false_in = (index == 0) & (ground_truth == 1); 
    true_out = (index == 1) & (ground_truth == 1); 
    false_out = (index == 1) & (ground_truth == 0); 
    
    
    c(true_in,:) = repmat(colors(2,:),sum(true_in==1),1);
    c(false_in,:) = repmat(colors(4,:),sum(false_in==1),1);
    c(true_out,:) = repmat(colors(1,:),sum(true_out==1),1);
    c(false_out,:) = repmat(colors(3,:),sum(false_out==1),1);
             
    scatter(clusterData(:,1), clusterData(:,2),150,c,'filled','MarkerFaceAlpha',0.7)
    
    title(['Embryo ', num2str(metaData.embryoID(embryoNumber))])
        
    xlabel min
    ylabel var
    
    grid on
    
    set(gca,'FontSize',12)
    
    embryoName = ['embryo_', sprintf('%03d', metaData.embryoID(embryoNumber))];
    
    savePath = fullfile(metaData.resultsFolder, embryoName);
    
    saveFigure(fig,saveName,savePath)
        
end


function  doAnalysis(results)

    stage = results.stage==3.0;
     
    ground_truth = repmat(results.ground_truth,1,4);
    index = results.classification;
           
    true_in = (index == 0) & (ground_truth == 0); 
    false_in = (index == 0) & (ground_truth == 1); 
    true_out = (index == 1) & (ground_truth == 1); 
    false_out = (index == 1) & (ground_truth == 0); 
    
    
    IN = grpstats((ground_truth == 0),results.embryo,'sum');
    OUT = grpstats((ground_truth == 1),results.embryo,'sum');
    
    % True In rate---------------------------------------------------------
    T_IN = grpstats(true_in,[results.embryo, stage],'sum');      
    TINR = T_IN./IN;    
    runStats(TINR,'true_in_rate',results.saveaPath)
    
    % In  miss rate--------------------------------------------------------
    F_OUT = grpstats(false_out,results.embryo,'sum');    
    FINR = F_OUT./IN;   
    runStats(FINR,'in_miss_rate',results.saveaPath)
    
    % True Out rate--------------------------------------------------------
    T_OUT = grpstats(true_out,results.embryo,'sum');     
    TOUTR = T_OUT./OUT;   
    runStats(TOUTR,'true_out_rate',results.saveaPath)
        
    % Out miss rate--------------------------------------------------------
    F_IN= grpstats(false_in,results.embryo,'sum');     
    FOUTR = F_IN./OUT;   
    runStats(FOUTR,'out_miss_rate',results.saveaPath)
           
end


function runStats(rates, saveName,savePath)

    colors = [1.0 0.5 0.4 ; % Salmon#
              0.4 0.6 0.9 ; % Cornflower Blue                          
              0.6 1.0 0.6 ; % Mint Green
              0.2 0.3 0.3]; %Charcoal
                    
    switch saveName
        case 'true_in_rate'
            Title = 'True Inside Rate';
            c = colors(2,:);
             figNum =2;
        case 'in_miss_rate'
            Title = 'Inside Miss Rate';
            c = colors(3,:);
            figNum =3;
        case 'true_out_rate'
             Title = 'True Outside Rate';
             c = colors(1,:);
             figNum =1;
        case 'out_miss_rate'
            Title = 'Outside Miss Rate';
            c = colors(4,:);
            figNum =4;
    end
          
    fig = figure(figNum);
    clf
    hold on
    boxplot(rates,'labels',{'MINS','Ellipsoidal','Convex Hull','insideOutside'},...
        'LabelOrientation','horizontal','Colors', c,'Symbol','k.');
    
    ylabel 'Rate'
    ylim([0,1.01])   
    title(Title)
    grid on
    
    set(gca, 'FontSize', 12,'YTick',0:0.25:1)
    
    ax = get(gca);
    
    for i = 1:numel(ax.Children.Children)
        ax.Children.Children(i).LineWidth = 2.5;
    end
    
    boxes = findobj(ax.Children.Children,'tag','Box');
    
    for i = 1:numel(boxes)
                
        x = boxes(i).XData;
        y = boxes(i).YData;
        
        h = fill(x,y,c,'FaceAlpha',0.5);
        
        uistack(h,'bottom')
        
       jitter(rates(:,i),i,c)
           
    end         
    
    shg
    
    [p,~,stats]=kruskalwallis(rates,[],'off');    
    multcompare(stats,'Display', 'off')    
    
    set(gca,'FontName','Arial','FontSize',15)
    set(gcf,'Position', [10 10 300 500]);
    saveStats(p, stats,saveName,savePath)
    
    saveFigure(fig,saveName,savePath)

%     pause(.1)
%     close(fig)
    
end


function jitter(rate,i,c)

    scatter(i+.25*(rand(size(rate))-0.5),rate,[],c,'filled','MarkerEdgeColor','k')

end


function saveFigure(fig,saveName,savePath)

    fullSaveName = fullfile(savePath, saveName);
    
    savefig(fig,fullSaveName)
    saveas(fig,fullSaveName,'svg')
    saveas(fig,fullSaveName,'png')

end


function saveStats(p,stats,saveName,savePath)
    
    saveName = strcat(saveName,'_stats');
    
end




















