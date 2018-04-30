function param_mat = genParamMat(pct_set,filt_dist_set,bail_thresh_set,...
    penalty_set,max_complex_set,noise_thresh_set,user_min_segs_set)

if nargin<8
    user_min_segs_set = [1];
end

param_mat = [];
 for pct_it = 1:length(pct_set)
     for filt_it = 1:length(filt_dist_set)
         for bt_it = 1:length(bail_thresh_set)
             for pen_it = 1:length(penalty_set)
                 for mc_it = 1:length(max_complex_set)
                     for nt_it = 1:length(noise_thresh_set)
                         for um_it = 1:length(user_min_segs_set)
                             param_mat = [param_mat; ...
                                 pct_set(pct_it) filt_dist_set(filt_it) ...
                                 bail_thresh_set(bt_it) penalty_set(pen_it) ...
                                 max_complex_set(mc_it) noise_thresh_set(nt_it) ...
                                 user_min_segs_set(um_it)];
                         end
                     end
                 end
             end
         end
     end
 end
end
