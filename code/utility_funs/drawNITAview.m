function drawNITAview(x_draw,y_draw,...
    tb,vi_type,...
    value_limits,doy_limits,date_limits,...
    draw_objid,info_col)

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

if nargin==7
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
all_OBJECTIDs = unique(tb.OBJECTID); 

if draw_objid == 9999 
    unique_OBJECTIDs = all_OBJECTIDs;
else 
    unique_OBJECTIDs = intersect(unique(draw_objid),all_OBJECTIDs);
end

vi_all = table2array(tb(:,{vi_type}));
for i=1:length(unique_OBJECTIDs)
    object_id = unique_OBJECTIDs(i);
    obj_idx = find(tb.OBJECTID==object_id); 
    
    if nargin==7
        temp = table2array(tb(:,{info_col}));
        plot_info = char(unique(temp(obj_idx)));
    else
        plot_info = '';
    end 
    
    
    im_date = date_num(obj_idx);
    doy = doy_all(obj_idx);
    vi = vi_all(obj_idx);
    
    [im_date,vi,doy] = filterLimits(im_date,vi,doy,value_limits,date_limits,doy_limits);
    
    handdraw_x = x_draw{i};
    handdraw_y = y_draw{i};
    
    figure, hold on
         scatter(im_date, vi, 100, doy,'.')
         plot(handdraw_x,handdraw_y,'red')
         colorbar;
         axis([min(im_date) max(im_date) value_limits])
         title([num2str(object_id) ' ' plot_info])
  
    pause(0.5);
end

end
