function [rmse_mean,rmse_med] = drawNITAcmp(tb,vi_type,...
    value_limits,doy_limits,date_limits,compute_mask,...
    param_mat,draw_objid,x_draw,y_draw)

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

[im_date_all,doy_all] = cellfun(@landsatImgDate,tb.system_index);
date_num = genDisDate(im_date_all,doy_all);
%sort table by date
table_sort = sortrows([array2table(date_num) tb], 1);
date_num = table_sort{:,1}; %converts back to double
tb = table_sort(:,2:end);

vi_all = table2array(tb(:,{vi_type}));

all_OBJECTIDs = unique(tb.OBJECTID); 
if draw_objid == 9999 
    unique_OBJECTIDs = all_OBJECTIDs;
else 
    unique_OBJECTIDs = intersect(unique(draw_objid),all_OBJECTIDs);
end

for param_it = 1:size(param_mat,2)
    param_line = param_mat{param_it};
    parfor i = 1:length(unique_OBJECTIDs)
        object_id = unique_OBJECTIDs(i);
        im_date = date_num(tb.OBJECTID==object_id);
        doy = doy_all(tb.OBJECTID==object_id);
        vi = vi_all(tb.OBJECTID==object_id);
        
        results_cell = nita_px(vi,im_date,doy,...
            value_limits,doy_limits,date_limits,param_line{1},param_line{2},...
            param_line{3},param_line{4},param_line{5},param_line{6},param_line{7},...
            compute_mask,param_line{8});
        
        x_nita{i} = results_cell{2};
        y_nita{i} = results_cell{3};
    
        date_vec_overlap_start = max(x_draw{i}(1),x_nita{i}(1));
        date_vec_overlap_end = min(x_draw{i}(end),x_nita{i}(end));
        
        draw_interp = interp1(x_draw{i},y_draw{i},date_vec_overlap_start:200:date_vec_overlap_end);
        nita_interp = interp1(x_nita{i},y_nita{i},date_vec_overlap_start:200:date_vec_overlap_end);

        % assess rmse of fit
        sq_error = (draw_interp-nita_interp).^2;
        rmse(i) = sqrt(mean(sq_error));  
    end
    rmse_mean(param_it) = mean(rmse);
    rmse_med(param_it) = median(rmse);
    
    if mod(param_it,10)==0
        disp(param_it)
    end
end

end






