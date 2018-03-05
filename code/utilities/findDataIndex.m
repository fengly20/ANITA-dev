% findDataIndex
% input: 
% return: non_nan_idx - a n*1 list taht contains indeces for good data
% points.
function non_nan_idx = findDataIndex(doy_limits, y, doy)
if isempty(doy_limits)
    non_nan_idx = find(not(isnan(y)));
else
    doy_filt = zeros(length(doy), 1);
    doy_limits_size = size(doy_limits); 
    for i = 1 : doy_limits_size(1) 
        doy_days = doy_limits(i, 1) : doy_limits(i, 2);
        doy_filt_temp = ismember(doy, doy_days);
        doy_filt = or(doy_filt, doy_filt_temp);
    end
    non_nan_idx = find(not(isnan(y))& doy_filt);
end  
end