
dataFolder = fullfile(pwd,'\figures\figure3\');

files = {'accuracy_insideOutside_raw_stats',...
         'accuracy_insideOutside_stats',...
         'accuracy_ellipsoidalMethod_stats',...
         'accuracy_convexHullMethod_stats'};

stat_types = {'sensitivity','specificity','precision'};
     
n_files = numel(files);
n_stats = numel(stat_types);

n_ball_radii = numel(r_ball);

summaryStats =  zeros(n_ball_radii,4,2);
summaryStats(:,1,1) = r_ball;
summaryStats(:,1,2) = r_ball;

for i = 1:n_files

    load(fullfile(dataFolder ,[files{i},'.mat']));
    
    for j = 1:n_stats
               
        [mu, sigma] = getSummaryStats( stats.(stat_types{j}));
        
        summaryStats(:,j+1,1) = mu;
        summaryStats(:,j+1,2) = sigma;
        
    end
        
    saveName = fullfile(dataFolder,[files{i},'.txt']);
    writeSummaryStats(summaryStats,saveName)
    
end


function [mu, sigma] = getSummaryStats(x)

    mu = mean(x,2);
    sigma =  std(x,[],2);
    
end


function writeSummaryStats(summaryStats,saveName)
    
    fileID = fopen(saveName,'w');
   
    str_format = '%15s %15s %15s %15s %15s %15s %15s %15s\n';
    str = {'r_Ball','mean_OUT','mean_IN','mean_PRECISION',...
           'r_Ball', 'std_OUT', 'std_IN', 'std_PRECISION'};
    
    fprintf(fileID,str_format,str{:});
   
    str_format = '%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f\n';
    
    for i = 1:size(summaryStats,1)   

        x = [summaryStats(i,:,1),summaryStats(i,:,2)];
        
        fprintf(fileID,str_format, x);

    end

   fclose(fileID);
   
end