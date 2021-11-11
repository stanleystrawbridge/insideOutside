function [stats] = calculateStats(idx, ground_truth)
    
    % Ground Truth counts 
    number_outside = sum(ground_truth == 1);
    number_inside = sum(ground_truth == 0);
    
    % Correct classification counts
    true_outside = sum((ground_truth == 1) & (idx == 1));    
    true_inside = sum((ground_truth == 0) & (idx == 0));    
    
    % Incorrect classification counts
    false_outside = sum((ground_truth == 0) & (idx == 1));    
    
    
    sensitivity = true_outside/number_outside; 
    specificity = true_inside/number_inside;
    precision = true_outside/(true_outside + false_outside);
    
    stats = [sensitivity, specificity, precision];
    
end
