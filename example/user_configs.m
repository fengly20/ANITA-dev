%% 
% VI parameters 
user_vi = 'ndvi';
vi_value_range = [-1 1];

%%
% NITA parameters 
%doy_limits = [1 366];
value_limits = vi_value_range;
doy_limits = [1 365]; %can enter in multiple DOY ranges to include 
date_limits = [-9999 9999];
bail_thresh = 1.2; 
noise_thresh = 1; 
penalty = 1.5; 
filt_dist = 7; 
pct = 90; 
max_complex = 10;
min_complex = 2;
compute_mask = 1;
filter_opt = 'movcv';

%%
% NITA metrics parameters
vi_change_thres = -0.3;
run_thres = 1000; % in days

%%
% NITA parameters for single pixel examination 
value_limits_dg = vi_value_range;
doy_limits_dg = [1 100; 260 365]; %can enter in multiple DOY ranges to include 
date_limits_dg = [1995000 2010000];
bail_thresh_dg = 1.2; 
noise_thresh_dg = 1; 
penalty_dg = 1; 
filt_dist_dg = 5; 
pct_dg = 50; 
max_complex_dg = 7;
min_complex_dg = 1;
compute_mask_dg = 1;
filter_opt_dg = 'movcv';
