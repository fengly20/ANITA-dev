%% Run the image stack 

%%
cd('D:/Feng/NITA-dev/control_file')
addpath(genpath('./'))
% add path for NITA 
addpath(genpath('../'))

%% 
%[vi_stack,vi_georeference] = geotiffread('./input/ak_tanana_NBR_stack_GEE.tif');
%dates_tb = readtable('./input/ak_tanana_dates_GEE.csv');
%create im_date and doy vectors from landsat system index column
%[im_date_all,doy_all] = cellfun(@landsatImgDate,dates_tb.system_index);

load('../input/vi_stack.mat')
load('../input/im_date_all.mat')
load('../input/doy_all.mat')

%distribute date across 1000
date_num = genDisDate(im_date_all,doy_all);
%%
% image_size_xy = (size(vi_stack,1)*size(vi_stack,2));
% for i = 1:size(vi_stack,3)
%     good_px_count = sum(sum(vi_stack(:,:,i)>-1)); %This value ("0") only works for NDVI!! But...
%     %a lot of clouds are not getting masked so not(isnan(pixel)) does not work well.
%     good_px_percent(i) = good_px_count/image_size_xy;
% end
% 
% %plot of good pixel percentage by image date
% figure, histogram(good_px_percent)
% 
% %find some good years to look at throughout your study's time range
% num_chunks = 6; %divide time into n chunks
% total_dist_days = max(date_num)-min(date_num);
% time_chunks = [min(date_num):total_dist_days/num_chunks:max(date_num)];
% 
% for i = 1:num_chunks
%    %find best image in each time chunk
%    image_chunk_idx = find(date_num>=time_chunks(i) & date_num<time_chunks(i+1));
%    best_image_per_chunk_idx = max(find(good_px_percent(image_chunk_idx)==...
%        max(good_px_percent(image_chunk_idx))))
%    best_images(i) = image_chunk_idx(best_image_per_chunk_idx);
% end
%    
% %display num_chunks figures
% figure('position',[200 200 1000 800]), hold on
% subplot1(2,3,'Gap',[0.02 0.02])
% for i = 1:num_chunks
%     subplot1(i)
%     %vi_contrast stretch
%     im_stretched = vi_stretch(vi_stack(:,:,best_images(i)));
%     imagesc(flipud(im_stretched)), axis equal;
%     title(num2str(date_num(best_images(i))))
% end
% 
% %% 
% %create a subset if necessary. In this workshop, running less than 100 px by 100 px is good.
% %Use the figure generated below to choose coordinates (do NOT use the figure from last section)
% im_stretched = vi_stretch(vi_stack(:,:,best_images(6)));
% figure, imagesc(im_stretched); axis equal;
% %xmin = 10; %specify upper left x
% %ymin = 109; %specify upper left y
% %xmax = 100; %lower right x
% %ymax = 199; %lower right y
% %vi_sub = vi_stack(ymin:ymax, xmin:xmax,:);
% 
% %get rid of images with less than 1% good px
% keep_images_dates_idx = find(good_px_percent>0.02);
% vi_sub = vi_stack(:,:,keep_images_dates_idx);
% date_num = date_num(keep_images_dates_idx);
% doy_all = doy_all(keep_images_dates_idx);

%%
vi_sub = vi_stack;

clear results_cells metrics_cell
results_cells = cell(size(vi_sub,1),size(vi_sub,2));
metrics_cells = cell(size(vi_sub,1),size(vi_sub,2));

parpool(5)

%%
user_configs
tic
for i = 1:size(vi_sub,1)
    parfor j = 1:size(vi_sub,2)
        vi_run = single(squeeze(vi_sub(i,j,:)));
        im_date_run = date_num;
        doy_run = doy_all;
        results_cell_run = nita_px(vi_run,im_date_run,doy_run,...
            value_limits,doy_limits,date_limits,bail_thresh,noise_thresh,...
            penalty,filt_dist,pct,max_complex,min_complex,...
            compute_mask,filter_opt);
        results_cells{i,j} = results_cell_run; 
    end
    i
end
toc

%%
%Unclear if parallel processing is beneficial here. Close call.
user_configs
tic
for i = 1:size(vi_sub,1)
    for j = 1:size(vi_sub,2)
        results_cell_run = results_cells{i,j};
        metrics_cell_run = computeMetrics(results_cell_run,vi_change_thres,run_thres,0.002739);
        metrics_cells{i,j} = metrics_cell_run; 
    end
    i
end
toc

%%
distdate_mat = distDateMat(metrics_cells,'middle');
dist_duration_mat = distDurationMat(metrics_cells); 
dist_mag_mat = distMagMat(metrics_cells); 
noise_mat = noiseMat(results_cells);
linerr_mat = linearErrorMat(results_cells);
norm_linerr_mat = errorNoiseRatioMat(results_cells);
complexity_mat = complexityMat(results_cells);
dist_slope_mat = distSlopeMat(metrics_cells);
post_dist_mat = postDistMat(metrics_cells);
total_value_change_mat = valueChangeMat(-9999,9999,'diff',metrics_cells);

duration = 2;
start_date = date2distdate(20000101);
end_date = duration*1000+start_date;
time_value_change_mat = valueChangeMat(start_date,end_date,'diff',metrics_cells,'progress_bar');

recovery2_mat = recoveryMat(2,'diff',metrics_cells,'progress_bar');
recovery10_mat = recoveryMat(10,'diff',metrics_cells,'progress_bar');
recovery2cmp_mat = recoverCmpMat(2,metrics_cells,'progress_bar');
[head_mat,tail_mat] = headTailMat(metrics_cells);

dispMat(distdate_mat);
dispMat(dist_duration_mat);
dispMat(dist_mag_mat);
dispMat(noise_mat);
dispMat(linerr_mat);
dispMat(norm_linerr_mat);
dispMat(complexity_mat);
dispMat(dist_slope_mat);
dispMat(post_dist_mat);
dispMat(total_value_change_mat);
dispMat(time_value_change_mat);
dispMat(recovery2_mat);
dispMat(recovery10_mat);
dispMat(recovery2cmp_mat);
dispMat(head_mat);
dispMat(tail_mat);

norm_rise = normRiseMat(results_cells); 
dispMat(norm_rise);
date = date2distdate(20000101);
flac_mat = flactuationMat(date,results_cells,metrics_cells,0.5);
dispMat(flac_mat);

%% 2.8 Reexamine certain pixels in metric images using updated NITA parameters
xi = 47;
yi = 73;

user_configs

%rerun NITA with diagnostic parameters (NOTE: This will not change your parameters that you set above)
vi_run = single(squeeze(vi_sub(yi,xi,:)));
im_date_run = date_num;
doy_run = doy_all;
results_cell_run = nita_px(vi_run,im_date_run,doy_run,...
    value_limits_dg,doy_limits_dg,date_limits_dg,bail_thresh_dg,noise_thresh_dg,...
    penalty_dg,filt_dist_dg,pct_dg,max_complex_dg,min_complex_dg,...
    compute_mask_dg,filter_opt_dg);

metrics_cell_run = computeMetrics(results_cell_run,vi_change_thres,run_thres,0.5);
%clf;
figure
[noise,fit_count] = viewNITA(vi_run,im_date_run,results_cell_run,doy_run,'fitvi','on')

%%
%http://geotiff.maptools.org/spec/geotiff6.html
if sum(size(vi_stack) ~= size(vi_sub))>0
    vi_georeference = geoRefSubset(vi_georeference,xmin,xmax,ymin,ymax);
end

proj_code = 32607;

geotiffwrite('./output/complexity.tif',complexity_mat,vi_georeference,'CoordRefSysCode',proj_code);
geotiffwrite('./output/distyear.tif',distyear_mat,vi_georeference,'CoordRefSysCode',proj_code);
geotiffwrite('./output/valuechange.tif',valuechange_mat,vi_georeference,'CoordRefSysCode',proj_code);
geotiffwrite('./output/datethres.tif',datethres_mat,vi_georeference,'CoordRefSysCode',proj_code);
geotiffwrite('./output/recover2.tif',recover2_mat,vi_georeference,'CoordRefSysCode',proj_code);
geotiffwrite('./output/recover4.tif',recover4_mat,vi_georeference,'CoordRefSysCode',proj_code);

geotiffwrite('./output/tanana_recov2.tif',recovery2_mat,vi_georeference,...
    'CoordRefSysCode',proj_code);

geotiffwrite('./output/tanana_val_change_full.tif',total_value_change_mat,vi_georeference,...
    'CoordRefSysCode',proj_code);

geotiffwrite('./output/tanana_dist_magnitude.tif',dist_mag_mat,vi_georeference,...
    'CoordRefSysCode',proj_code);