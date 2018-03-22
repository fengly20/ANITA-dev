function [metrics_cell] = computeMtrics(results_cell,vi_change_thres,run_thres,recovery, vi_thres)
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
% recov_4yr - 4 years recovery value after first disturbance  
% recov_2yr - 2 years recovery value after first disturbance  
%%
% ---
% 1. extract information from results_cell
  knots = results_cell{2};
  coeffs = results_cell{3};
  rises = results_cell{7};
  runs = results_cell{8};
  runs_in_days = results_cell{9};

try    
    change_percent = rises./coeffs(1:end-1);
    slopes = rises./runs*1000; %1000 is one dist_date year
 
  % 2. disturbance detection 
    dist_idx = find(change_percent<vi_change_thres & runs_in_days<=run_thres);
  % 2.a case of disturbance exists     
    if not(isempty(dist_idx))  
        num_dist = length(dist_idx); % this is an output 
        cum_mag_dist =sum(abs(rises(dist_idx))); % this is an output; cumulative magnitude disturbance
        % min_idx = find(coeffs==min(coeffs),1);
        
      % 2.a.1 metrics for first disturbance             
        dist_date_before = knots(dist_idx(1)); % this is an output
        dist_date_nadir = knots(dist_idx(1)+1); % this is an output
        coeff_nadir = coeffs(dist_idx(1)+1); % this is an output
      %this ONLY works with distributed dates!
        dist_duration = dist_date_nadir - dist_date_before; % this is an output
        dist_slope = slopes(dist_idx(1)); % this is an output
        %dist_mag = rises(dist_idx(1));
      
      % 2.a.2 first "recovery" (or just segment after first disturbance)
        if length(slopes)>dist_idx(1) % this is the case that the disturbance is not the last segament
            post_dist_slp = slopes(dist_idx(1)+1); % this is an output
            post_dist_mag = rises(dist_idx(1)+1); % this is an output
        else % this is the case that the disturbance is the last segament
            post_dist_slp = NaN; % this is an output
            post_dist_mag = NaN; % this is an output
        end
      
      % 2.a.3 recovery stuff (optional)
        if nargin > 3
          % get some useful values 
            study_year_first = knots(1); 
            study_year_last = knots(end);   
            study_year_first_str = num2str(study_year_first);
            study_year_last_str = num2str(study_year_last);
            study_years = str2double(study_year_first_str(1:4)):str2double(study_year_last_str(1:4)); 
          
          % interpolation 
            study_years_1000 = study_years*1000;
            study_years_1000(end) = study_years_1000(end)+999;
            valid_years_idx = find(study_years_1000>=knots(1)-999 & study_years_1000<=knots(end)+999);
            coeff_interp = interp1(knots,coeffs,study_years_1000(valid_years_idx),'linear');
            coeffs_by_year(valid_years_idx) = coeff_interp;
            coeffs_by_year(1) =coeffs(1);
            coeffs_by_year(end) = coeffs(end);
        
          % find the recovery value   
            [~,temp_idx]= min(abs(dist_date_nadir+2000-study_years_1000));
            %approximate_2yrs_date = study_years_1000(temp_idx); 
            recov_2yr = coeffs_by_year(temp_idx);
            
            [~,temp_idx]= min(abs(dist_date_nadir+4000-study_years_1000));
            %approximate_4yrs_date = study_years_1000(temp_idx); 
            recov_4yr = coeffs_by_year(temp_idx);
            
          % find the vi_thres date
            [~,temp_idx]= min(abs(coeffs_by_year-vi_thres));
            if min(abs(coeffs_by_year-vi_thres))<(0.1*vi_thres)
                thres_year = floor(study_years_1000(temp_idx)/1000);
            else 
                thres_year = NaN;
            end
        else
            recov_2yr = NaN;
            recov_4yr = NaN;
            thres_year = NaN;
        end % end of if nargin > 3            

  % 2.b case that there's no disturbance   
    else
        num_dist = NaN;
        cum_mag_dist = NaN;
        dist_date_before =NaN;
        dist_date_nadir =NaN;
        dist_duration =NaN;
        dist_slope =NaN;
        coeff_nadir = NaN;
        post_dist_slp =NaN;
        post_dist_mag=NaN;

        recov_4yr = NaN;
        recov_2yr = NaN;
        
        thres_year = NaN;
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

    recov_4yr = NaN;
    recov_2yr = NaN;
    
    thres_year = NaN;
end % end of try 
                  
metrics_cell = {num_dist cum_mag_dist...
    dist_date_before dist_date_nadir...
    dist_duration dist_slope...
    coeff_nadir...
    post_dist_slp post_dist_mag...
    recov_2yr recov_4yr...
    thres_year};     
end % end of the function 
       