function [x_draw, y_draw] = drawNita(tb,vi_type,value_limits,doy_limits,date_limits,draw_objid,info_col)

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
    
    figure, hold on
         scatter(im_date, vi, 100, doy,'.')
         colorbar;
         axis([min(im_date) max(im_date) value_limits])
         title([num2str(object_id) ' ' plot_info])

     [x_draw{i}, y_draw{i}] = ginput();
     
    close
    
    % perform a correction on hand-draw trajectory 
    % possible corrections are:
    % (1) correct double click (non unique data pairs)
    % (2) correct backwards x values (assuming y values are reasonable)

    % (1)
    [~,~,tp_uni_idx] = unique([x_draw{i} y_draw{i}],'row','stable');
    x_draw{i} = x_draw{i}(tp_uni_idx);
    y_draw{i} = y_draw{i}(tp_uni_idx);

    % (2)
    for j=2:length(x_draw{i})
        tp_bw_flag = (x_draw{i}(j) - x_draw{i}(j-1)) <= 0;
        if tp_bw_flag
            x_draw{i}(j) = x_draw{i}(j-1)+3;
        end
    end
    
end

end
