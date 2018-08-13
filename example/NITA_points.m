%% Run NITA on points data 

%%
% set path 
cd('/Users/leyang/Documents/AU-projects/NITA-dev/example')
% add path 
addpath(genpath('./'))
addpath(genpath('../code'))

%%
pts_table = readtable('./input/test_pts.csv');

%%   
[im_date_all,doy_all] = cellfun(@landsatImgDate,pts_table.system_index);

%distribute date across 1000
date_num = genDisDate(im_date_all,doy_all);

% add all iinfo into the table 
pts_table = [array2table(date_num) array2table(doy_all) pts_table];

%% 1.4 Run NITA on a single point and view the plot
user_configs

object_id = 1;

pt_objectid = pts_table(pts_table.OBJECTID==object_id,:);
pt_objectid = sortrows(pt_objectid,1);

im_date = pt_objectid.date_num;
doy = pt_objectid.doy_all;
%nbr = pt_objectid.nbr;
ndvi = pt_objectid.ndvi;

if strcmp(user_vi,'ndvi')
    vi = ndvi;
elseif strcmp(user_vi,'nbr')
    vi = nbr;
end

%profile on
results_cell = nita_px(vi,im_date,doy,...
    value_limits,doy_limits,date_limits,bail_thresh,noise_thresh,...
    penalty,filt_dist,pct,max_complex,min_complex,...
    compute_mask,filter_opt);
%profile viewer
%[bail_cut,fit_count] = viewNITA(vi,im_date,results_cell,doy,'allvi','on')
figure
[bail_cut,fit_count] = viewNITA(vi,im_date,results_cell,doy,'fitvi','on')

% %% 1.6 Run NITA on all points (max=25) and visualize together 
% % Now we can run NITA on all of the provided points and view the plot all together. 
% %%
% site_ids = sort(unique(pts_table.OBJECTID))';
% 
% subplot_ncol = min(ceil(sqrt(length(site_ids))),5);
% subplot_nrow = min(ceil(length(site_ids)/subplot_ncol),5);
% bail_cuts = zeros(1,length(site_ids));
% 
% figure, subplot1(subplot_nrow,subplot_ncol,'Gap',[0.02 0.02])
% 
% for i = 1:min((subplot_ncol*subplot_nrow),length(site_ids))
%     
%     object_id = site_ids(i);
%     
%     pt_objectid = pts_table(pts_table.OBJECTID==object_id,:);
%     pt_objectid = sortrows(pt_objectid,1);
% 
%     im_date = pt_objectid.date_num;
%     doy = pt_objectid.doy_all;
%     %nbr = pt_objectid.nbr;
%     ndvi = pt_objectid.ndvi;
%     
%     if strcmp(user_vi,'ndvi')
%         vi = ndvi;
%     elseif strcmp(user_vi,'nbr')
%         vi = nbr;
%     end
%     
%     user_configs
%     
%     results_cell = nita_px(vi,im_date,doy,...
%         value_limits,doy_limits,date_limits,bail_thresh,noise_thresh,...
%         penalty,filt_dist,pct,max_complex,min_complex,...
%         compute_mask,filter_opt);
%         
%     subplot1(i);   
%     [bail_cut,fit_count] = viewNITA(vi,im_date,results_cell,doy,'fitvi','off');
%     title([num2str(object_id) '  ' num2str(bail_cut)])
%     bail_cuts(1,i) = bail_cut;
% end
% 
% %% 
% max(bail_cuts)
% min(bail_cuts)
% mean(bail_cuts)
