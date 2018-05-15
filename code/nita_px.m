function results_cell = nita_px(px,date_vec,doy_vec,...
    value_limits,doy_limits,date_limits,bail_thresh,noise_thresh,...
    penalty,filt_dist,pct,max_complex,min_complex,...
    compute_mask,filter_opt)
%% Input arguments: 
% Data: 
%     'px'
%     'date_vec'
%     'doy_vec'
% Constraints: 
%     'value_limits'
%     'doy_limits'
%     'date_limits'
%     'bail_thresh'
%     'noise_thresh'
% Numerical args: 
%     'penalty'
%     'filt_dist'
%     'pct'
%     'max_complex'
%     'min_complex'
% Switches: 
%     'compute_mask'
% Options: 
%     'filter_opt'

%% Documentation 
%anita code purpose:
%using a time series of spectral information (e.g., NDVI, NBR) generate a
%single-line or piecewise fit by adding breakpoints at x,y locations of
%change. This algorithm is "insensitive" to noise in that it: 1) filters
%out noise pixel-dates based on a user threshold, 2) uses noise-adaptive
%filtering to determine the optimal location of break points. That is, a
%breakpoint will be more likely placed in a region of change with low noise
%than high noise because we have more confidence that the change is real;
%3) user enters pct which can "float" the fit on top of data that are noisy
%in the downward direction (e.g., NDVI). This can also be used to fit
%annual peak greenness, avoiding phenological variation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%image_line: submit your image line by line. This facilitates
%parallelization.

%date_vec: vector of image dates. Dates should be in the form of... CHECK
%ON THIS

%penalty: penalty parameter for Bayesian Information Criterion (BIC). A
%higher penalty will lead to fewer segments in the piecewise fit. Range is
%(generally) 1-10 where 1 allows more segments and 10 will almost certainly
%give you a single-line fit.

%bail_thresh: User set parameter determing whether or not to run the full
%nita code or if a linear fit is adequate or if the data are so noisy that
%we can't hope to do better than a linear fit. Bail_thresh is compared to
%mae_lin/noise. If there is high linear error and low noise this ratio will
%be high, indicating a disturbance, thus the code should run. Bail thresh
%of 1 or below is conservative in that it will mostly run the full code.

%max_complex: set the maximum number of segments you'll allow in your
%piecewise fit. 10 is a good upper limit for most applications. This does
%NOT mean that you will end up with 10 segments. Some will be removed based
%on BIC.

%filt_dist: Set the size of the search window (i.e., number of adjacent
%image dates) over which to search for breaks. Higher filt_dist values
%(e.g., 7 or 9) will result in fewer sharp breaks. Higher is better in very
%noisy data but 3 is frequently a good starting place.

%pct: "percentile" -- this is the "float" parameter. How do you want your
%fit to float on the data? If you're modeling NDVI change, you may want pct
%of 75 or 90 to float your fit at the 75th or 90th percentile because the
%least contaminated NDVI values are generally the highest. With non-noisy
%NBR, a pct of 50 (right down the middle) may be fine.

%doy: day-of-year vector from 001 to 365 same length as im_date

%doy_limits: optional parameter where user can specify the doy range to
%include in the fit. For instance, in Alaska, I only use doy in summer
%months to establish the fit (e.g., [170 230])

%noise thresh: User decides how much noise is too much. Noise here is
%defined as the forward finite difference in the spectral values.
%Basically, how much are adjacent values jumping around. This setting is
%dependent on the range of values in the input dataset. For NBR data
%ranging from 0 to around 6500, I used a value of 2000 for noise_thresh.

%diag_plots: do you want the diagnostic plots to pop up when you run this
%code? The answer (if you're calling this code as a function) is probably:
%NO. If you run a full image with diag_plots == 1, your computer will
%crash. The right answer is almost always "0".

%%
% ---
% 0. check the inputs

% check the data inputs 

% if input image line is not double, it must be converted
  if ~isa(px,'double')
      px = double(px);
  end
  
% if dates come in as 1xn, need to transpose to nx1 to match the output of squeeze
  if size(px,2)>size(px,1)
      px = px';
  end
  
  if size(date_vec,2)>size(date_vec,1)
      date_vec = date_vec';
  end
  
  if size(doy_vec,2)>size(doy_vec,1)
      doy_vec = doy_vec';
  end
  
% remove possible duplicated im_date
% for example it's possible to get two different values for the same
% distributed date (due to image overlaping), in such case the first value
% is kept (no pericular reason to decide which one is discarded). 
  [~, unq_idx] = unique(date_vec);
  date_vec = date_vec(unq_idx);
  px = px(unq_idx);
  doy_vec = doy_vec(unq_idx);

%%  
  try
%%      
% --- 
% 0.5 prepare x and y 
      x = date_vec;
      y = px;
    
    % apply value_limits, doy_limits and date_limits 
      [x,y,doy_vec] = filterLimits(x,y,doy_vec,value_limits,date_limits,doy_limits);
            
    %noise calc (in spectral index units)
      noise = median(abs(diff(y)));
            
    %filter by abs(diff) to get rid of bad noise per user threshold
      diff_holder = diff(y);
      good_idx = find(abs(diff_holder) <= noise_thresh)+1;
      x = x(good_idx);
      y = y(good_idx);
      x_len = length(x);

% ---
% 0.6 sanity check of x and y 
    % check for adequate date 
      if x_len <= (filt_dist*2) 
          error('Not enough data pairs!'); 
      end
%%
% ---
% 1. single line fit 
    %set starting coeffs (first date and last date and first VI and last VI)
      first_coeff = prctile(y(1:filt_dist),pct);
      last_coeff = prctile(y(end-filt_dist:end),pct);

      knot_set =  [x(1);x(end)]; 
      coeff_set = [first_coeff;last_coeff];
      coeff_indices = [1;x_len];
      
      pts = [x y];

    % calculate ortho error for the single line fit 
      dist_init = calDistance(knot_set,coeff_set, pts);
      mae_lin = calMae(dist_init);

    % diagnostic plot for linear fit section
      %figure, hold on 
      %plot(x,y,'.')
      %plot(knot_set, coeff_set)
          
%%      
% ---
% 2. NITA
      if mae_lin/noise > bail_thresh && compute_mask==1 %determine whether to run full code
%%       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %In here, is the nita build phase, where breakpoint
        %locations are selected and added up to max_complex
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
        % set starting conditions for for-loop
          mae_ortho(1) = mae_lin;
          
        % this will run until max_complex or until there are no more 
        % aviable breakpoints can be added .
          for i = 2:max_complex 
              % ortho error using the current knot set
                clear dist
                dist = calDistance(knot_set,coeff_set,pts);
                [cand_idx,coeff] = findCandidate(dist,filt_dist,pct,y,coeff_indices,filter_opt);
         
                if cand_idx == -999
                    break
                end
              
                [knot_set,coeff_set,coeff_indices] = updateknotcoeffSet(knot_set,coeff_set,coeff_indices,x,cand_idx,coeff);
                dist_new = calDistance(knot_set,coeff_set, pts);
                mae_ortho(i) = calMae(dist_new);
          end % end of for i=2:max_complex  
 
          complexity_count = length(knot_set)-1;
          
        % grab final error
          mae_final = mae_ortho(complexity_count);
        
          %figure, hold on
          %plot(x,y,'r.')
          %plot(knot_set,coeff_set,'o')
          %axis([min(x) max(x) -5000 7000])
  
%%          
% ---
% 3. BIC removal process 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %now take knots away iteratively
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
        % *_max saved as copies 
          knots_max = knot_set;
          coeffs_max = coeff_set;
        
        % keep_* will go into for loop and shrink every time   
          keep_knots = knot_set;
          keep_coeffs = coeff_set;
          %mae_ortho(1) = mae_final

          yinterp1 = interp1(knots_max,coeffs_max,x,'linear');
          y_pos_idx = (y-yinterp1)>0;  
          
          for i=1:complexity_count-(min_complex-1)
            % loop through knots, removing each and checking which
            % raises MAE the least compared to weighted data
              keep_idx{i} = genKeepIdx(keep_knots,keep_coeffs,pts,pct,y_pos_idx);
   
            % need to reccalc ortho_err BEFORE the if statement
            % because the bic will now be the point of entry
              clear dist
              dist = calDistance(keep_knots,keep_coeffs,pts);
              ortho_err = min(dist,[],2);
              mae_ortho_holder(i) = calMae(dist);
              
            % here is the reweighting of the error based on pct. If
            % you chose pct = 75 then places where your fit
            % underestimates are more important to fix.
              ortho_err(y_pos_idx) = ortho_err(y_pos_idx)*pct;
              ortho_err(~y_pos_idx) = ortho_err(~y_pos_idx)*(100-pct);
    
              bic_remove(i) = calBIC(ortho_err,keep_knots,penalty);  

              if i==1
                  knot_storage{i} = knots_max;
                  coeff_storage{i} = coeffs_max;
              else                       
                  knot_storage{i} = keep_knots(keep_idx{i});
                  coeff_storage{i} = keep_coeffs(keep_idx{i});

            % reduce the number of knots and coeffs each
            % iteration
                  keep_coeffs = keep_coeffs(keep_idx{i});
                  keep_knots = keep_knots(keep_idx{i});
              end 
          end
              
          bic_idx = find(bic_remove==min(bic_remove));
          keep_coeffs = coeff_storage{bic_idx};
          keep_knots = knot_storage{bic_idx};
          mae_final = mae_ortho_holder(bic_idx);    
             
          %figure, hold on
          %plot(x,y,'r.')
          %plot(keep_knots,keep_coeffs,'o')
          %plot(keep_knots, keep_coeffs)
          %axis([min(x) max(x) -5000 7000])
          
        %outputs assuming the NITA code ran
          complexity = length(keep_knots)-1;
          final_knots = keep_knots;
          final_coeffs = keep_coeffs;
        %error (mae_ortho)
          mae_final_ortho = mae_final;
          mae_linear = mae_lin;
          noise_out = noise;
      else
          complexity = 1;
          final_knots = [min(x);max(x)];
          final_coeffs = [first_coeff; last_coeff];
          mae_final_ortho = mae_lin;
          mae_linear = mae_lin;
          noise_out = noise;
      end % end of  if mae_lin/noise <= bail_thresh && compute_mask == 1
      rises = diff(final_coeffs);
      runs = diff(final_knots);     
      runs_days = runs/1000*365;
  catch %error_report 
      complexity = -999;
      final_knots = -999;
      final_coeffs = -999;
      mae_final_ortho = -999;
      mae_linear = -999;
      noise_out = -999;
      rises = -999;
      runs = -999;
      runs_days = -999;
      pts = -999;
  end % end of try-catch

%results output
  results_cell = {complexity final_knots final_coeffs mae_linear mae_final_ortho noise_out rises runs runs_days pts};

% print out the error message
%  if exist(error_report,'var') == 1
%      disp(error_report.message)
%  end
  
end %end of the function
        