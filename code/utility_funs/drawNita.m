function [x_draw, y_draw] = drawNita(tb,vi_type,vi_value_range,doy_limits,max_draw,info_col)

% ---
%1. check the table
column_names = lower(tb.Properties.VariableNames);
if ismember('objectid',column_names)~= 1
    error('No OBJECTID column in the table');
end

if ismember(vi_type,column_names)~= 1
    error(['No ' vi_type ' column in the table']);
end

if ismember('system_index',column_names)~= 1
    error('No system_index column in the table');
end

if nargin==6
    if ismember(info_col,column_names)~= 1
        error(['No ' info_col ' column in the table']);
    end 
end

[im_date_all,doy_all] = cellfun(@landsatImgDate,tb.system_index);
date_num = genDisDate(im_date_all,doy_all);
%sort table by date
table_sort = sortrows([array2table(date_num) tb], 1);
date_num = table_sort{:,1}; %converts back to double
tb = table_sort(:,2:end);

% ---
%2. loop through each point 
unique_OBJECTIDs = unique(tb.OBJECTID);
draw_count = min(length(unique_OBJECTIDs),max_draw);
vi_all = table2array(tb(:,{vi_type}));
for i=1:draw_count
    object_id = unique_OBJECTIDs(i);
    obj_idx = find(tb.OBJECTID==object_id); 
    
    if nargin==6
        temp = table2array(tb(:,{info_col}));
        plot_info = char(unique(temp(obj_idx)));
    else
        plot_info = '';
    end 
    
    
    im_date = date_num(obj_idx);
    doy = doy_all(obj_idx);
    vi = vi_all(obj_idx);
    
    non_nan_idx = findDataIndex(doy_limits, vi, doy);
    vi = vi(non_nan_idx);
    im_date = im_date(non_nan_idx);
    doy = doy(non_nan_idx);
    
    vi_idx = find(vi>vi_value_range(1) & vi<=vi_value_range(end));
    vi = vi(vi_idx);
    im_date = im_date(vi_idx);
    doy = doy(vi_idx);
    
    figure, hold on
         scatter(im_date, vi, 100, doy,'.')
         colorbar;
         axis([min(im_date) max(im_date) vi_value_range])
         title([num2str(object_id) ' ' plot_info])

     [x_draw{i}, y_draw{i}] = ginput();
     
    close
end

end







