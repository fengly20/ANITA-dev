function [results_cell] = nita_px(px, date_vec, penalty,...
    bail_thresh, max_complex, filt_dist, pct, doy, doy_limits,...
    noise_thresh, diag_plots)

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


% ---
% 0. check the input px, date_vec

%if input image line is not double, it must be converted
  if ~isa(px,'double')
      px = double(px);
  end

%if dates come in as 1xn, need to transpose to nx1 to match the output
%of squeeze
  if size(date_vec,2)>size(date_vec,1)
      date_vec = date_vec';
  end
    
%define the cell outputs in advance of loop. Final_knots stores date
%values and final_coeffs stores spectral index values of breakpoints.
%  final_knots = cell(size(px,2),1);
%  final_coeffs = cell(size(px,2),1);

%loop through all samples of line (1 x numSamp x numDates)
  
  try
      
% --- 
% 0.5 prepare x and y 
    %set x and y
      x = date_vec;
      y = px;
           
    %screen or do not screen for seasonal limits using doy_limit
    %parameter
      non_nan_idx = findDataIndex(doy_limits, y, doy);
      x = double(x(non_nan_idx));
      y = double(y(non_nan_idx));
            
    %noise calc (in spectral index units)
      noise = median(abs(diff(y)));
            
    %filter by abs(diff) to get rid of bad noise per user threshold
      diff_holder = diff(y);
      good_idx = find(abs(diff_holder) <= noise_thresh)+1;
      x = x(good_idx);
      y = y(good_idx);
      x_len = length(x);

% ---
% 1. single line fit 
    %set starting coeffs (first date and last date and first SI
    %value and last SI value)
      first_coeff = prctile(y(1:filt_dist),pct);
      last_coeff = prctile(y(end-filt_dist:end),pct);

    %use the line established by the first_coeff and last_coeff above
    %use orthogonal distance from that line to find location of
    %first candidate breakpoint
    %ortho error
      knot_set =  [min(x);max(x)];
      coeff_set = [first_coeff;last_coeff];
      coeff_indices = [1;x_len];
      
      pts = [x y];

      dist_init = calDistance(knot_set,coeff_set, pts);
      mae_lin = calMae(dist_init);

    %set knot location and establish tentative knot set (comprised
    %currently of start date, first candidate breakpoint, and end
    %date)
      %knot_loc = x(cand_idx);
      %[knot_set,coeff_set,coeff_indices] = updateknotcoeffSet(knot_set,coeff_set,coeff_indices,knot_loc,cand_idx,coeff);
           
    %diagnostic plot for linear fit section
      if diag_plots==1
          figure, hold on 
          plot(x,y,'.')
          plot([min(x) max(x)],[first_coeff last_coeff])
          plot(knot_set, coeff_set)
          %plot(x,mov_cv*1,'-g')
      end
      
% ---
% 2. NITA
      if mae_lin/noise > bail_thresh %determine whether to run full code
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %In here, is the nita build phase, where breakpoint
        %locations are selected and added up to max_complex
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
        %set starting conditions for while loop
          complexity_count = length(knot_set)-1;
          mae_ortho(complexity_count) = mae_lin;
        %this will run until max_complex or until there are no more
        %viable breakpoints to add.
          while complexity_count~=max_complex 
              complexity_count = complexity_count + 1;
            %ortho error using the current knot set
              clear dist
              dist = calDistance(knot_set,coeff_set, pts);
              [cand_idx,coeff,search_series] = findCandidate(dist,filt_dist,pct,y);
                    
              while sum(ismember(coeff_indices,cand_idx))~=0   
                  [cand_idx,coeff] = findNextCandidate(coeff_indices,search_series,filt_dist,pct,y);        
              end
         
              if cand_idx == -999
                  break
              end
              
              [knot_set,coeff_set,coeff_indices] = updateknotcoeffSet(knot_set,coeff_set,coeff_indices,x,cand_idx,coeff);
              dist_new = calDistance(knot_set,coeff_set, pts);
              mae_ortho(complexity_count) = calMae(dist_new);
              
              %diagnostic plotting
%             figure, hold on
%             plot(x,y,'.')
%             plot(knot_set_hold,coeff_set,'or')
%             axis([min(x) max(x) -2 2])
          end % end of while complexity_count~=max_complex 
 
        %grab final error
          mae_final = mae_ortho(complexity_count-1);
        %working on roll back
          %keep_knots = knots_prev;
          %keep_coeffs = coeffs_prev;
          keep_knots = knot_set;
          keep_coeffs = coeff_set;
             
          if diag_plots == 1
              figure, hold on
              plot(x,y,'r.')
              plot(keep_knots,keep_coeffs,'o')
              line(keep_knots, keep_coeffs)
%             axis([min(x) max(x) -5000 7000])
          end % end of if diag_plots == 1
          
          
% ---
% 3. now take knots away iteratively
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %now take knots away iteratively
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          clear mae_ortho; clear bic_remove; clear aic_remove; clear num_segs;
          clear coeff_storage; clear knot_storage;
        
          knots_max = keep_knots;
          coeffs_max = keep_coeffs;
          %mae_ortho(1) = mae_final;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Do while there are more than 2 knots being evaluated. That
        %is, check the BIC for num knots = max_complex
        %incrementally down to num knots = 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          pos=1;
          [~, unq_idx] = unique(keep_knots);
          yinterp1 = interp1(knots_max(unq_idx),coeffs_max(unq_idx),x,'linear');
          y_pos_idx = (y-yinterp1)>0;  
          while pos ~= complexity_count

            %loop through knots, removing each and checking which
            %raises MAE the least compared to weighted data
              keep_idx{pos} = genKeepIdx(keep_knots,keep_coeffs,pts,pct,y_pos_idx);
   
            %need to reccalc ortho_err BEFORE the if statement
            %because the bic will now be the point of entry
              clear dist
              dist = calDistance(keep_knots,keep_coeffs,pts);
              ortho_err = min(dist,[],2);
            %here is the reweighting of the error based on pct. If
            %you chose pct = 75 then places where your fit
            %underestimates are more important to fix.
              ortho_err(y_pos_idx) = ortho_err(y_pos_idx)*pct;
              ortho_err(~y_pos_idx) = ortho_err(~y_pos_idx)*(100-pct);
              mae_ortho_holder(pos) = mean(min(dist,[],2));
            
            [num_segs(pos),bic_remove(pos)] = calBIC(ortho_err,keep_knots,penalty);  

            if pos==1
                knot_storage{1} = knots_max;
                coeff_storage{1} = coeffs_max;
            else                       
                knot_storage{pos+1} = keep_knots(keep_idx{pos});
                coeff_storage{pos+1} = keep_coeffs(keep_idx{pos});

              %reduce the number of knots and coeffs each
              %iteration
                keep_coeffs = keep_coeffs(keep_idx{pos});
                keep_knots = keep_knots(keep_idx{pos});
            end % end of if pos==1
              pos = pos + 1;
          end % end of while pos ~= complexity_count
   
          
          bic_idx = find(bic_remove==min(bic_remove));
          keep_coeffs = coeff_storage{bic_idx};
          keep_knots = knot_storage{bic_idx};
          mae_final = mae_ortho_holder(bic_idx);    
             
          if diag_plots == 1
              figure, hold on
              plot(x,y,'r.')
              plot(keep_knots,keep_coeffs,'o')
              ine(keep_knots, keep_coeffs)
%             axis([min(x) max(x) -5000 7000])
          end % end of if diag_plots == 1
          
        %outputs assuming the NITA code ran
          complexity = length(keep_knots)-1;
          final_knots = keep_knots;
          final_coeffs = keep_coeffs;
        %error (mae_ortho)
          mae_final_ortho = mae_final;
          mae_linear = mae_lin;
          noise_out = noise;
      end % end of if mae_lin/noise > bail_thresh
    %only populate for linear model if not already populated by
    %NITA
      if mae_lin/noise <= bail_thresh
          complexity = 1;
          final_knots = [min(x);max(x)];
          final_coeffs = [first_coeff; last_coeff];
          %error (mae)
          mae_final_ortho = mae_lin;
          mae_linear = mae_lin;
          noise_out = noise;
      end % end of  if mae_lin/noise <= bail_thresh
      rises = diff(final_coeffs);
      runs = diff(final_knots);     
      runs_days = runs/1000*365;
  catch
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

%"progress bar" (can reinstate if using multicore)
%       ff('line: %u\n',line_num)
%       fprintf('avg pixel time for this line: %u\n',mean(toc_i))
 
end %end of function
        

