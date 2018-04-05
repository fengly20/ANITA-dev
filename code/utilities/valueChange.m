function out = valueChange(start_date,end_date,option,metrics_cell)

% start_date: -9999, or distributed date number 
% duration: 9999, or distributed date number
% option: 'diff', 'change_percent', 'value'
try 
    interp_pts = metrics_cell{10};

    x_dates = interp_pts(:,1);
    y_vis = interp_pts(:,2);
    
    if start_date >= x_dates(1) && start_date <= x_dates(end)
        [~,x_date_idx] = min(abs(x_dates-start_date));
        x_date_idx = x_date_idx(1);
        x_start_date = x_dates(x_date_idx);
    end

    if end_date > x_start_date && end_date <= x_dates(end)
        [~,x_date_idx] = min(abs(x_dates-end_date));
        x_date_idx = x_date_idx(end);
        x_end_date = x_dates(x_date_idx);
    end 
    
  % grab value and calculate change 
    y_start_date = y_vis(x_dates == x_start_date);
    y_end_date = y_vis(x_dates == x_end_date);
    
    if strcmp(option,'diff')
        out = y_end_date - y_start_date;
    elseif strcmp(option,'change_percent')
        out = (y_end_date - y_start_date)/y_start_date*100;
    elseif strcmp(option,'value')
        out = y_end_date;
    else
        error('non-valid option provided!')
    end

catch
    out = NaN;
end
end

    
    