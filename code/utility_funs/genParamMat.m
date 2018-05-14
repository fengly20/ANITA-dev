function param_mat = genParamMat(bail_thresh_set,noise_thresh_set,penalty_set,filt_dist_set,pct_set,max_complex_set,min_complex_set,filter_opt_set)

param_mat = {};
count = 0;
for bt = 1:length(bail_thresh_set)
    for nt = 1:length(noise_thresh_set)
        for pen = 1:length(penalty_set)
            for fd = 1:length(filt_dist_set)
                for pct = 1:length(pct_set)
                    for maxc = 1:length(max_complex_set)
                        for minc = 1:length(min_complex_set)
                            for fo = 1:length(filter_opt_set)
                                 temp_cell = ...
                                     {bail_thresh_set(bt) ...
                                     noise_thresh_set(nt) ...
                                     penalty_set(pen) ...
                                     filt_dist_set(fd) ...
                                     pct_set(pct) ...
                                     max_complex_set(maxc) ...
                                     min_complex_set(minc) ...
                                     filter_opt_set(fo)};
                                 count = count + 1;
                                 param_mat{count} = temp_cell;
                            end
                        end
                    end
                end
            end
        end
    end
end

end
