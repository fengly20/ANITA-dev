function [x_draw, y_draw] = drawNitaAmend(tb,...
    vi_type,value_limits,doy_limits,date_limits,...
    draw_objid,redo_objid,x_draw,y_draw,info_col)

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

if nargin==8
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
%2. find redo id location 
if draw_objid == 9999 
    draw_objid = unique(tb.OBJECTID);
end

loc_idx = find(ismember(draw_objid,redo_objid));

% ---
%3. redo the hand draw 
if isempty(loc_idx)
    disp('No OBJECTID provided to re draw!')
else
    vi_all = table2array(tb(:,{vi_type}));
    for i=1:length(loc_idx)
        tp_locidx = loc_idx(i);
        object_id = draw_objid(tp_locidx);
        obj_idx = find(tb.OBJECTID==object_id); 
    
        if nargin==8
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

        [x_draw{tp_locidx}, y_draw{tp_locidx}] = ginput();
     
        close
    
        % perform a correction on hand-draw trajectory 
        % possible corrections are:
        % (1) correct double click (non unique data pairs)
        % (2) correct backwards x values (assuming y values are reasonable)

        % (1)
        [~,~,tp_uni_idx] = unique([x_draw{tp_locidx} y_draw{tp_locidx}],'row','stable');
        x_draw{tp_locidx} = x_draw{tp_locidx}(tp_uni_idx);
        y_draw{tp_locidx} = y_draw{tp_locidx}(tp_uni_idx);

        % (2)
        for j=2:length(x_draw{tp_locidx})
            tp_bw_flag = (x_draw{tp_locidx}(j) - x_draw{tp_locidx}(j-1)) <= 0;
            if tp_bw_flag
                x_draw{tp_locidx}(j) = x_draw{tp_locidx}(j-1)+3;
            end
        end
    
    end
end

end
