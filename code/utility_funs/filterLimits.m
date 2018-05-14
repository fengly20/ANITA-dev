function [x,y,doy_vec] = filterLimits(x,y,doy_vec,value_limits,date_limits,doy_limits)

% the filter is in four fold: 
% (1) no x values out of dates_limits
% (2) no x values out of doy_limits
% (3) no nan y values 
% (4) no y values out of value_range  

% (1)
if sum(date_limits) == 0 % date_limits = [-9999 9999]
    x = x;
    y = y;
    doy_vec = doy_vec; 
else % all other cases 
    if date_limits(1)==-9999 && date_limits(2)~=9999
        start_date = x(1);
        end_date = date_limits(2);
    elseif date_limits(1)~=-9999 && date_limits(2)==9999
        start_date = date_limits(1);
        end_date = x(end);
    else
        start_date = date_limits(1);
        end_date = date_limits(end);
    end
    
    date_idx = find(x>=start_date & x<=end_date);
    x = x(date_idx);
    y = y(date_idx);
    doy_vec = doy_vec(date_idx);
end

% (2)
doy_idx = zeros(length(doy_vec),1);
doy_limits_size = size(doy_limits); 
for i = 1 : doy_limits_size(1) 
    doy_days = doy_limits(i, 1) : doy_limits(i, 2);
    doy_idx_temp = ismember(doy_vec, doy_days);
    doy_idx = or(doy_idx, doy_idx_temp);
end
x = x(doy_idx);
y = y(doy_idx);
doy_vec = doy_vec(doy_idx);

% (3)
non_nan_idx = find(not(isnan(y)));
x = x(non_nan_idx);
y = y(non_nan_idx);
doy_vec = doy_vec(non_nan_idx);

% (4)
value_idx = find(y>value_limits(1) & y<=value_limits(end));
x = x(value_idx);
y = y(value_idx);
doy_vec = doy_vec(value_idx);

end
