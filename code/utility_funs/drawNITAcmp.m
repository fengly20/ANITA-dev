function [rmse_mean,rmse_med] = drawNITAcmp(tb,vi_type,vi_value_range,doy_limits,param_mat,max_draw,x_draw,y_draw,filter_opt)

if exist('filter_opt','var')==0
    filter_opt = 'movcv';
end

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

unique_OBJECTIDs = unique(tb.OBJECTID);
vi_all = table2array(tb(:,{vi_type}));
for params_it = 1:size(param_mat,1)
    param_line = param_mat(params_it,:);
    parfor i = 1:min(max_draw,length(unique(tb.OBJECTID)))
        
        object_id = unique_OBJECTIDs(i);
        obj_idx = find(tb.OBJECTID==object_id); 
    
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

        results_cell = nita_px(vi,im_date,param_line(4),param_line(3),...
            param_line(5),param_line(2),...
            param_line(1),doy,doy_limits,param_line(6),0,param_line(7),filter_opt);
        
        %[bail_cut,fit_count] = viewNITA(vi,im_date,results_cell,doy,'fitvi','on')
        
        x_nita{i} = results_cell{2};
        y_nita{i} = results_cell{3};

        date_vec_overlap_start = max(x_draw{i}(1),x_nita{i}(1));
        date_vec_overlap_end = min(x_draw{i}(end),x_nita{i}(end));
        draw_interp = interp1(x_draw{i}, y_draw{i},[date_vec_overlap_start:200:date_vec_overlap_end]);
        nita_interp = interp1(x_nita{i}, y_nita{i},[date_vec_overlap_start:200:date_vec_overlap_end]);

        %assess rmse of fit
        sq_error = (draw_interp-nita_interp).^2;
        rmse(i) = sqrt(mean(sq_error));
     end
     rmse_mean(params_it) = mean(rmse);
     rmse_med(params_it) = median(rmse);
     if mod(params_it,10)==10
         params_it
     end
end

end






