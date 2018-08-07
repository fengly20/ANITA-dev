%% Run NITA on points data 

%%
% set path 
cd('D:/Feng/NITA-dev/control_file')
% add path 
addpath(genpath('./'))
addpath(genpath('../code'))

%%
pts_table = readtable('../input/pts2.csv');

%% 
[im_date_all,doy_all] = cellfun(@landsatImgDate,pts_table.system_index);

%distribute date across 1000
date_num = genDisDate(im_date_all,doy_all);

% add all iinfo into the table 
pts_table = [array2table(date_num) array2table(doy_all) pts_table];

%% 1.4 Run NITA on a single point and view the plot
user_configs

object_id = 3;

pt_objectid = pts_table(pts_table.OBJECTID==object_id,:);
pt_objectid = sortrows(pt_objectid,1);

im_date = pt_objectid.date_num;
doy = pt_objectid.doy_all;
%nbr = pt_objectid.nbr;
ndvi = pt_objectid.ndvi;

if strcmp(user_vi,'ndvi')
    vi = ndvi;
end


user_configs
results_cell = nita_px(vi,im_date,doy,...
    value_limits,doy_limits,date_limits,bail_thresh,noise_thresh,...
    penalty,filt_dist,pct,max_complex,min_complex,...
    compute_mask,filter_opt);
%[bail_cut,fit_count] = viewNITA(vi,im_date,results_cell,doy,'allvi','on')
figure
[bail_cut,fit_count] = viewNITA(vi,im_date,results_cell,doy,'fitvi','on')

%%

clear x_draw y_draw

draw_objid = [1 3 5 7 9];

user_configs
[x_draw,y_draw] = drawNita(pts_table,user_vi,value_limits,doy_limits,date_limits,draw_objid);

%save('./control_file/x_draw','x_draw')
%save('./control_file/y_draw','y_draw')
%% 
if not(exist('x_draw','var'))
    load('./control_file/x_draw','x_draw')
end

if not(exist('y_draw','var'))
    load('./control_file/y_draw','y_draw')
end

drawNITAview(x_draw,y_draw,...
    pts_table,user_vi,...
    value_limits,doy_limits,date_limits,...
    draw_objid)

%% Only run this section if any of the hand draw needs to be re-do
redo_objid = [3 7];

[x_draw,y_draw] = drawNitaAmend(pts_table,...
    user_vi,value_limits,doy_limits,date_limits,...
    draw_objid,redo_objid,x_draw,y_draw);

%%
%gather stats for line comparisons

bail_thresh_set = [1.3 1.5 1.7 1.9 2.1];
noise_thresh_set = [1];
penalty_set = [0.5 1 2 3 5];
filt_dist_set = [1 3 5 7];
pct_set = [50 70 80 90];
max_complex_set = [5 7 10 15];
min_complex_set = [1];
filter_opt_set = ["movcv","zscore","movmean","gaussian","lowess","loess","rlowess","rloess","sgolay"];

%assemble parameter mat
param_mat = genParamMat(bail_thresh_set,noise_thresh_set,penalty_set,filt_dist_set,pct_set,max_complex_set,min_complex_set,filter_opt_set);

%%
if isempty(gcp('nocreate'))
    parpool(4)
end
    
clear rmse_mean rmse_med pct95_mean
user_configs
[rmse_mean,rmse_med,pct95_mean] = drawNITAcmp(pts_table,user_vi,...
    value_limits,doy_limits,date_limits,compute_mask,...
    param_mat,draw_objid,x_draw,y_draw);

%% 
% How should we filter the parameter combinations? 
low_err_idx = find(pct95_mean<=min(pct95_mean));
selected_param_mat = {param_mat{low_err_idx}};

%% 1.6 Run NITA on all points (max=25) and visualize together 
% Now we can run NITA on all of the provided points and view the plot all together. 
%%
% choose a param line
param_line = selected_param_mat{1};

site_ids = sort(unique(pts_table.OBJECTID))';

subplot_ncol = min(ceil(sqrt(length(site_ids))),5);
subplot_nrow = min(ceil(length(site_ids)/subplot_ncol),5);
bail_cuts = zeros(1,length(site_ids));

figure, subplot1(subplot_nrow,subplot_ncol,'Gap',[0.02 0.02])

for i = 1:min((subplot_ncol*subplot_nrow),length(site_ids))
    
    object_id = site_ids(i);
    pt_objectid = pts_table(pts_table.OBJECTID==object_id,:);
    pt_objectid = sortrows(pt_objectid,1);

    im_date = pt_objectid.date_num;
    doy = pt_objectid.doy_all;
    %nbr = pt_objectid.nbr;
    ndvi = pt_objectid.ndvi;
    
    user_configs
    results_cell = nita_px(vi,im_date,doy,...
            value_limits,doy_limits,date_limits,param_line{1},param_line{2},...
            param_line{3},param_line{4},param_line{5},param_line{6},param_line{7},...
            compute_mask,param_line{8});
        
    subplot1(i);   
    [bail_cut,fit_count] = viewNITA(vi,im_date,results_cell,doy,'fitvi','off');
    title([num2str(object_id) '  ' num2str(bail_cut)])
    bail_cuts(1,i) = bail_cut;
end

%% 
format short g
fprintf('bail_thresh = %i \n',param_line{1});...
    fprintf('\n');...
    fprintf('noise_thresh = %i \n',param_line{2});...
    fprintf('\n');...
    fprintf('penalty = %i \n',param_line{3});...
    fprintf('\n');...    
    fprintf('filt_dist = %i \n',param_line{4});...
    fprintf('\n');...
    fprintf('pct = %i \n',param_line{5});...
    fprintf('\n');... 
    fprintf('max_complex = %i \n',param_line{6});...
    fprintf('\n');...
    fprintf('min_complex = %i \n',param_line{7});...
    fprintf('\n');...
    fprintf('filter_opt = %s \n',param_line{8});...
    fprintf('\n');
