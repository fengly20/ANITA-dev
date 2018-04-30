function dist_date = calDistDate(metrics_cell,option)

% three valid values for option:
% (1) 'begining'
% (2) 'middle'
% (3) 'end'
 
dist_date_before = metrics_cell{3};
dist_date_nadir = metrics_cell{4};

if strcmp(option,'beginning')
    dist_date = dist_date_before;
elseif strcmp(option,'middle')
    dist_date = (dist_date_before+dist_date_nadir)/2;
elseif strcmp(option,'end')
    dist_date = dist_date_nadir;
else 
    error('non-valid option is given! use ''begining'', ''middle'', ''end''.')
end
end