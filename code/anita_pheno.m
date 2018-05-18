function [results_cell] = anita_pheno(image_px, date_vec, penalty,...
    bail_thresh, max_complex, filt_dist, pct, doy, date_range,...
    noise_thresh, diag_plots)

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
%image_px: submit your image line by line. This facilitates
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
%noisy data but 3 is freqently a good starting place.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %if input image line is not double, it must be converted
    if ~isa(image_px,'double')
        image_px = double(image_px);
    end

    %if dates come in as 1xn, need to transpose to nx1 to match the output
    %of squeeze
    if size(date_vec,2)>size(date_vec,1)
        date_vec = date_vec';
    end
    
    %define the cell outputs in advance of loop. Final_knots stores date
    %values and final_coeffs stores spectral index values of breakpoints.
    final_knots = cell(size(image_px,2),1);
    final_coeffs = cell(size(image_px,2),1);

        %loop through all samples of line (1 x numSamp x numDates)
 
        try
            
            clear mad; clear mov_cv; clear slm; clear mae_ortho;
            clear coeff_indices;
            
            %set x and y
            x = doy;
            y=squeeze(image_px);
           
            %find only those points in the date range specified
            if isempty(date_range)
                non_nan_idx = find(not(isnan(y)));
            else
                non_nan_idx = find(not(isnan(y)) & date_vec>date_range(1) &...
                    date_vec<date_range(2));
            end
            x = double(x(non_nan_idx));
            y = double(y(non_nan_idx));
                                    
            %sort data by DOY
            xy = [x y];
            xy_sort = sortrows(xy,1);
            x = xy_sort(:,1); y=xy_sort(:,2);
%             year_sort = sortrows([year_vec_sub xy],2);
%             year_sort = year_sort(:,1);
            
            %noise calc (in spectral index units)
            noise = median(abs(diff(y)));
            
            %filter by abs(diff) to get rid of bad noise per user threshold
            diff_holder = diff(y);
            good_idx = find(diff_holder>-noise_thresh)+1;
            x = x(good_idx);
            y = y(good_idx);

            %set starting coeffs (first date and last date and first SI
            %value and last SI value)
            first_coeff = prctile(y(1:filt_dist),pct);
            last_coeff = prctile(y(end-filt_dist:end),pct);

            %use the line established by the first_coeff and last_coeff above
            %use orthogonal distance from that line to find location of
            %first candidate breakpoint
            %ortho error
             pts = [x y];
             edge = [min(x) first_coeff...
             max(x) last_coeff];
             dist_init = distancePointEdge(pts, edge);

            %assess noiseness in local area (filt_dist) using moving inverse cv
            mov_mean = movmean(dist_init,filt_dist);
            mov_std = movstd(dist_init, filt_dist);
            mov_cv = mov_mean./mov_std;
           
            %mean absolute orthogonal error from simple linear fit. Used to
            %determine whether the piecewise portion of the code should run
            %or if a single-segment fit is adequate.
            mae_lin = mean(dist_init);                   

            %Seek location of first candidate breakpoint based on maximum
            %of mov_cv
            cand_idx = find(abs(mov_cv)==max(abs(mov_cv)),1);
           %find local filt_dist nhood of cand_idx
            cand_idx_filt = cand_idx-((filt_dist-1)/2):1:cand_idx+((filt_dist-1)/2);
            
            %if the candidate point is with filt_dist of beginning or end
            %of time series, ignore and search for another point
            if cand_idx > length(y) - filt_dist || cand_idx < filt_dist
                %find the highest regional maxima of mov_cv and select that
                %instead of the selection above
                pks_bin = imregionalmax(mov_cv);
                mov_cv_sort = unique(sort(mov_cv(pks_bin)));
                cand_idx = find(mov_cv==mov_cv_sort(end-1),1,'last');
                cand_idx_filt = cand_idx-((filt_dist-1)/2):1:cand_idx+((filt_dist-1)/2);
                coeff = prctile(y(cand_idx_filt),pct);
            else
                %if there is no issue, go ahead and set the coeff at the
                %chosen breakpoint based on pct and candidate's nhood
                coeff = prctile(y(cand_idx_filt),pct);
            end
            
            %set knot location and establish tentative knot set (comprised
            %currently of start date, first candidate breakpoint, and end
            %date)
            knot_loc = x(cand_idx);
            knot_set = [min(x); knot_loc; max(x)];
            
            %make corresponding first set of coeffs
            coeff_set = [first_coeff; coeff; last_coeff];
            %keep track of coeff indices so I can easily insert new ones
            coeff_indices = [1; cand_idx; length(y)];
           
            %diagnostic plot for linear fit section
            if diag_plots==1
                figure, hold on 
                    plot(x,y,'.')
                    plot([min(x) max(x)],[first_coeff last_coeff])
                    plot(x,mov_cv*100,'-g')
            end
%%
            if mae_lin/noise > bail_thresh %determine whether to run full code
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %In here, is the nita build phase, where breakpoint
                %locations are selected and added up to max_complex
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                %set starting conditions for while loop
                pos = 1;        
                not_done = true;
                %this holds the error quantity at each build step. Starts
                %with the linear error from above.
                mae_ortho(pos) = mae_lin;
                
                %this will run until max_complex or until there are no more
                %viable breakpoints to add.
                while not_done
   
                    %hold initial knots and coeffs for rollback
                    %also hold previous set of mad points for comparision
                    if pos>1
                        knots_prev = knot_set_hold;
                        mov_cv_prev = mov_cv;
                    else
                        knots_prev = [min(x) max(x)];
                        mov_cv_prev = mov_cv;
                    end

                    %hold the current set of knots for reversion. Necessary
                    %because knot_set will have new knot added
                     knot_set_hold = knot_set;
                    
                     %increment pos now because we've already completed one
                     %iteration (essentially) by establishing the
                     %single-segment fit
                     pos = pos+1;
                             
                    %ortho error using the current knot set
                     pts = [x y];
                    clear dist
                    for m = 1:length(knot_set)-1
                        edge = [knot_set(m) coeff_set(m)...
                            knot_set(m+1) coeff_set(m+1)];
                        dist(:,m) = distancePointEdge(pts, edge);
                    end
                    ortho_err = min(dist,[],2);
                    
                    %This doesn't actually matter during build, right?
                    mae_ortho(pos) = median(ortho_err);
                     
                     %diagnostic plotting
%                      figure, hold on
%                         plot(x,y,'.')
%                         plot(knot_set_hold,coeff_set,'or')
%                         %axis([min(x) max(x) -2 2])

                                  
                   %if we've reached max_complex then wrap things up. We
                   %will have overshot max_complex so, need to roll
                   %everything back by one iteration
                   if pos==max_complex 
                        %break out of the loop
                        not_done=false;
                        %grab final error
                        mae_final = mae_ortho(pos-1);

                        %working on roll back
                        keep_knots = knots_prev;
                        keep_coeffs = coeffs_prev;
                   else

                        %again, use inverse cv to separate signal from
                        %noise
                        mov_mean = movmean(ortho_err,filt_dist);
                        mov_std = movstd(ortho_err, filt_dist);
                        mov_cv = mov_mean./mov_std;
                        
                        %find new cand_idx based on max mov_cv
                        %i think this cv_change is unused
                        %cv_change = abs(mov_cv_prev)-abs(mov_cv);
                        cand_idx = find(mov_cv==max(mov_cv),1,'last');

                        %if this new candidate has NOT previously been
                        %selected, proceed, else, need to search for next
                        %most viable candidate.
                        if sum(ismember(coeff_indices,cand_idx))==0
                            %update knot set
                            knot_loc = x(cand_idx);  
                            knot_set = sort([knot_set; knot_loc]);
                        else
                            dupe_cand_idx = true;
                            counter = 1;
                            while dupe_cand_idx
                                pks_bin = imregionalmax(mov_cv);
                                mov_cv_sort = unique(sort(mov_cv(pks_bin)));
                                cand_idx = find(mov_cv==mov_cv_sort(end-counter),1,'last');
                                counter = counter +1;
                                if sum(ismember(coeff_indices,cand_idx))==0 ||...
                                        counter>length(mov_cv_sort)-1
                                    dupe_cand_idx = false;
                                end
                            end
                            knot_loc = x(cand_idx);  
                            knot_set = sort([knot_set; knot_loc]);
                        end
                   end %pos==max_complex if statement
                                       
                            
                    %assign coeff based pct y value of nearby points
                    cand_idx_filt = cand_idx-((filt_dist-1)/2):1:cand_idx+((filt_dist-1)/2);
                    if cand_idx > length(y) -filt_dist
                        coeff = prctile(y((cand_idx-filt_dist):cand_idx),pct);
                    elseif cand_idx < filt_dist
                        coeff = prctile(y(cand_idx:(cand_idx+filt_dist)),pct);
                    else
                        coeff = prctile(y(cand_idx_filt),pct);
                    end
                    %store previous before adding new
                    coeffs_prev = coeff_set;
                    
                    %insert new coeff in appropriate place in sequence
                    coeff_set_unsort = [coeff; coeff_set];
                    coeff_indices_unsort = [cand_idx; coeff_indices];
                    coeffs_and_indices_sorted = sortrows([coeff_set_unsort coeff_indices_unsort],2);
                    coeff_set = coeffs_and_indices_sorted(:,1);
                    coeff_indices = coeffs_and_indices_sorted(:,2);

                end %end of while loop
                
                 if diag_plots == 1
                    figure, hold on
                            plot(x,y,'r.')
                            plot(keep_knots,keep_coeffs,'o')
                            line(keep_knots, keep_coeffs)
                            axis([min(x) max(x) -5000 7000])
                end
%%                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %now take knots away iteratively
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                clear mae_ortho; clear bic_remove; clear aic_remove; clear num_segs;
                clear coeff_storage; clear knot_storage;
                knots_max=keep_knots;
                coeffs_max = keep_coeffs;
                mae_ortho(1) = mae_final;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Do while there are more than 2 knots being evaluated. That
                %is, check the BIC for num knots = max_complex
                %incrementally down to num knots = 2
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                pos=1;
                not_done = true;
                while not_done 
                    
                    %with all knots build an interplated representation of the fit 
                    %with a point at each image date. This is assumed to be
                    %the most "accurate" representation of the fit to the
                    %data and thus provides a baseline against which to
                    %assess the lowered accuracy of removing breakpoints.
                    if pos == 1
                        [~, unq_idx] = unique(keep_knots);
                        yinterp1 = interp1(keep_knots(unq_idx),keep_coeffs(unq_idx),x,'linear');
                    end
                    
                    %show the interpolated line in red
%                     if diag_plots==1
%                         figure, hold on
%                             plot(x,y,'k.')
%                             plot(x, yinterp1,'.r')
%                             plot(keep_knots,keep_coeffs,'o')
%                             %axis([min(x) max(x) -3.5 1.5])
%                     end


                   %loop through knots, removing each and checking which
                   %raises MAE the least compared to weighted data
                   clear mae_iter_ortho
                   if length(keep_knots)-2 > 0
                       
                       %isolate data points above interp fit. 
                       y_pos_idx = (y-yinterp1)>0;
                       
                       for j = 1:length(keep_knots)-2 %only inner knots

                           remove_idx = j+1;
                           include_idx = setdiff(1:length(keep_knots),remove_idx);

                            new_knot_set = keep_knots(include_idx);

                            coeff_set = keep_coeffs(include_idx);

                            %compute new mae_ortho
                             pts = [x y];
                            clear dist
                            for m = 1:length(new_knot_set)-1
                                edge = [new_knot_set(m) coeff_set(m)...
                                    new_knot_set(m+1) coeff_set(m+1)];
                                dist(:,m) = distancePointEdge(pts, edge);
                            end
                            ortho_err = min(dist,[],2);
                            ortho_err(y_pos_idx)=ortho_err(y_pos_idx)*pct;
                            ortho_err(~y_pos_idx) = ortho_err(~y_pos_idx)*(100-pct);
                            mae_iter_ortho(j+1) = mean(ortho_err);

                       end %end of iterative cv_change loop

                        %now find index for actual knot removal based on
                        %minimizing the mae_iter_ortho error
                        remove_idx = find(mae_iter_ortho==...
                            min(mae_iter_ortho(2:end)),1);
                        %keep all knots EXCEPT for the remove_idx
                        keep_idx{pos} = setdiff(1:length(keep_knots),remove_idx);
                   else
                       %this will only be first and last knots
                       keep_idx{pos} = [1 3];
                   end
                   
                    %need to reccalc ortho_err BEFORE the if statement
                    %because the bic will now be the point of entry
                    pts = [x y];
                    clear dist
                    for m = 1:length(keep_knots)-1
                        edge = [keep_knots(m) keep_coeffs(m)...
                            keep_knots(m+1) keep_coeffs(m+1)];
                        dist(:,m) = distancePointEdge(pts, edge);
                    end
                    ortho_err = min(dist,[],2);
                    %here is the reweighting of the error based on pct. If
                    %you chose pct = 75 then places where your fit
                    %underestimates are more important to fix.
                    ortho_err(y_pos_idx)=ortho_err(y_pos_idx)*pct;
                    ortho_err(~y_pos_idx) = ortho_err(~y_pos_idx)*(100-pct);
                    mae_ortho_holder(pos) = mean(min(dist,[],2));
                   
                   %BIC acccumulation
                    positive_idx = find(ortho_err>0); %in case a value is exactly 0
                    %this line is maybe the most computationally expensive
                    %single line in anita
                    dist_ortho_err = fitdist(ortho_err(positive_idx), 'lognormal');
                    loglik = -dist_ortho_err.negloglik;
                    num_segs(pos)=length(keep_knots)-1;

                    %user-defined penalty parameter
                    bic_remove(pos) = -2*loglik+penalty*num_segs(pos)*log(length(ortho_err));


                     %only ending loop when all added knots are removed
                    if pos == max_complex-1
                        not_done = false;
                       
                        %select best model based on min bic
                        bic_idx = find(bic_remove==min(bic_remove));                        
                       
                        keep_coeffs = coeff_storage{bic_idx};
                        keep_knots = knot_storage{bic_idx};
                        mae_final = mae_ortho_holder(bic_idx);
                    else
                        %hold the knots and coeffs from each iteration
                        %For the pos==1 run store both all the knots from
                        %the build-up phase AND the first set with
                        %one-removed (pos+1)
                        if pos==1
                            knot_storage{1} = knots_max;
                            coeff_storage{1} = coeffs_max;
                        end
                        
                        knot_storage{pos+1} = keep_knots(keep_idx{pos});
                        coeff_storage{pos+1} = keep_coeffs(keep_idx{pos});

                        %reduce the number of knots and coeffs each
                        %iteration
                        keep_coeffs = keep_coeffs(keep_idx{pos});
                        keep_knots = keep_knots(keep_idx{pos});

                    end

                    pos = pos + 1;

                end
               
                if diag_plots == 1
                    figure, hold on
                            plot(x,y,'r.')
                            plot(keep_knots,keep_coeffs,'o')
                            line(keep_knots, keep_coeffs)
                            axis([min(x) max(x) -5000 7000])
                end
%%
               
                %outputs assuming the NITA code ran
                    complexity = length(keep_knots)-1;
                    final_knots = keep_knots;
                    final_coeffs = keep_coeffs;
                    %error (mae_ortho)
                    mae_final_ortho = mae_final;
                    mae_linear = mae_lin;
                    noise_out = noise;
               

            else %else for bail thresh if

            %only populate for linear model if not already populated by
            %NITA         
                complexity = 1;
                final_knots = [x(1); x(end)];
                final_coeffs = [first_coeff; last_coeff];
                %error (mae)
                mae_final_ortho = mae_lin;
                mae_linear = mae_lin;
                noise_out = noise;
           
            
            end

                        
        catch
            complexity = -999;
            final_knots = -999;
            final_coeffs = -999;
            mae_final_ortho = -999;
            mae_linear = -999;
            noise_out = -999;
        end %end of try-catch
        

    
    %results output
     results_cell = {complexity final_knots final_coeffs...
         mae_linear mae_final_ortho noise_out};

      %"progress bar" (can reinstate if using multicore)
%       fprintf('line: %u\n',line_num)
%       fprintf('avg pixel time for this line: %u\n',mean(toc_i))
    
    
end %end of function
        

