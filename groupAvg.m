function group = groupAvg(data)
    % returns -1 if decay, 0 if sustained, 1 = growth
    t = length(data);
    data = data.^2;
    
    section_1 = mean(data(1:round(t/4)));
    section_2 = mean(data(round(t/4)+1:2*round(t/4)));
    section_3 = mean(data(2*round(t/4)+1:3*round(t/4)));
    section_4 = mean(data(3*round(t/4)+1:t));
    
    sec = [section_1 section_2 section_3 section_4]';
    diff = [-1 1 0 0;0 -1 1 0;0 0 -1 1]*sec;
    
    threshold = 100;
    diff_sum = sum(diff);
        
    group = diff;
    
    
end