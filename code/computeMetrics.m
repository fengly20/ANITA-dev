function [metrics_cell] = computeMetrics(results_cell,vi_change_thres,run_thres,time_step)
%%
% The function for now only cares about first disturbance in the time
% series, any scenarios assuming having multiple disturbances may end up
% with getting some metrics that do not make sense at all. 
% The list below contains all metrics calculated by this function: 
% num_dist - number of disturbances
% cum_mag_dist - cumulative magnitude of the disturbance(s)
% dist_date_before - the start date of first disturbance
% dist_date_nadir - the date that has the lowest value for that disturbance
% dist_duration - the duration of the first disturbance (time duration is also called 'run' in the function)
% dist_slope - the slope of the disturbance (which will be always a negative number)
% coeff_nadir - the lowest value of that disturbance
% post_dist_slp - slope of the next segament right after the first disturbance 
% post_dist_mag - magnitude of the next segament right after the first disturbance 
% 
%%
% ---
% 0. check number of input arguments 
if nargin ~= 4
    error('Not enough input arguments.\n ')
end

% ---
% 1. extract information from results_cell
knots = results_cell{2};
coeffs = results_cell{3};
rises = results_cell{7};
runs = results_cell{8};
runs_in_days = results_cell{9};

try    
    % 1.5 interpolating points 
    % get some useful values 
    study_year_first = knots(1); 
    study_year_last = knots(end);   
    study_year_first_str = num2str(study_year_first);
    study_year_last_str = num2str(study_year_last);
    study_years = str2double(study_year_first_str(1:4)):str2double(study_year_last_str(1:4)); 
          
    % interpolation 
    study_years_dis = [study_years(1):time_step:(study_years(end)+1)]*1000;
    valid_years_idx = study_years_dis>knots(1) & study_years_dis<knots(end);
    study_years_all = [knots(1) study_years_dis(valid_years_idx) knots(end)];
    coeff_interp = interp1(knots,coeffs,study_years_all,'linear');
    interp_pts = [study_years_all' coeff_interp'; knots(2:end-1) coeffs(2:end-1)];
    interp_pts_test = unique(interp_pts,'row'); % this is an output  
    
  % 2. disturbance detection
    change_percent = rises./abs(coeffs(1:end-1));
    slopes = rises./runs*1000; %1000 is one dist_date year
    dist_idx = find(change_percent<vi_change_thres & runs_in_days<=run_thres);
  % 2.a case of disturbance exists     
    if not(isempty(dist_idx))  
        num_dist = sum(diff(dist_idx)>1)+1; % this is an output
        % num_dist = length(dist_idx); % this is an output 
        cum_mag_dist =sum(abs(rises(dist_idx))); % this is an output; cumulative magnitude disturbance
         
      % 2.a.1 metrics for first disturbance             
        dist_date_before = knots(dist_idx(1)); % this is an output
        if length(dist_idx)>1 && (dist_idx(2)-dist_idx(1))==1
            dist_date_nadir = knots(dist_idx(2)+1); % this is an output
            coeff_before = coeffs(dist_idx(1)); % this is an output
            coeff_nadir = coeffs(dist_idx(2)+1); % this is an output
          %this ONLY works with distributed dates!
            dist_duration = dist_date_nadir-dist_date_before; % this is an output
            dist_mag = abs(coeff_nadir-coeff_before); % this is an output
            dist_slope = dist_mag/dist_duration; % this is an output
               
            if length(slopes)>dist_idx(2) % this is the case that the disturbance is not the last segament
                post_dist_slp = slopes(dist_idx(2)+1); % this is an output
                post_dist_mag = rises(dist_idx(2)+1); % this is an output
            else % this is the case that the disturbance is the last segament
                post_dist_slp = NaN; % this is an output
                post_dist_mag = NaN; % this is an output
            end
                    
        else % length(dist_idx)==1 or more disturbance event
            dist_date_nadir = knots(dist_idx(1)+1); % this is an output
            coeff_before = coeffs(dist_idx(1)); % this is an output
            coeff_nadir = coeffs(dist_idx(1)+1); % this is an output
          %this ONLY works with distributed dates!
            dist_duration = dist_date_nadir-dist_date_before; % this is an output
            dist_slope = slopes(dist_idx(1)); % this is an output
            dist_mag = rises(dist_idx(1)); % this is an output
      
          % 2.a.2 first "recovery" (or just segment after first disturbance)
            if length(slopes)>dist_idx(1) % this is the case that the disturbance is not the last segament
                post_dist_slp = slopes(dist_idx(1)+1); % this is an output
                post_dist_mag = rises(dist_idx(1)+1); % this is an output
            else % this is the case that the disturbance is the last segament
                post_dist_slp = NaN; % this is an output
                post_dist_mag = NaN; % this is an output
            end
        end     
        
  % 2.b case that there's no disturbance   
    else
        num_dist = NaN;
        cum_mag_dist = NaN;
        dist_date_before = NaN;
        dist_date_nadir = NaN;
        dist_duration = NaN;
        dist_slope = NaN;
        coeff_nadir = NaN;
        post_dist_slp = NaN;
        post_dist_mag = NaN;
        dist_mag = NaN;
        coeff_before = NaN;

    end % end of if not(isempty(dist_idx))  
    
catch %in the case of no data at all
    num_dist = NaN;
    cum_mag_dist = NaN;
    dist_date_before = NaN;
    dist_date_nadir = NaN;
    dist_duration = NaN;
    dist_slope = NaN;
    coeff_nadir = NaN;
    post_dist_slp = NaN;
    post_dist_mag = NaN;
    interp_pts = NaN;
    dist_mag = NaN;
    coeff_before = NaN;
    
end % end of try 
                  
metrics_cell = {num_dist cum_mag_dist...
    dist_date_before dist_date_nadir...
    dist_duration dist_slope...
    coeff_nadir...
    post_dist_slp post_dist_mag...
    interp_pts dist_mag coeff_before};     
end % end of the function 
       