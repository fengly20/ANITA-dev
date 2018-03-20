function [cand_idx,coeff] = findNextCandidate(coeff_indices,search_series_current,filt_dist,pct,y)
  
  search_series_n = search_series_current;
  search_series_n(coeff_indices) = [];
  
  
  if length(search_series_n) == (filt_dist-1)*2
      cand_idx = -999;
      coeff = -999;
  else
      cand_idx = find(search_series_current==max(search_series_n),1);
      
      not_done = true;
      search_series_less = search_series_n;
      cand_idx_less = find(search_series_n==max(search_series_n),1);
      while not_done
          if cand_idx > (length(y)-filt_dist) || cand_idx < (filt_dist+1)
              search_series_less(cand_idx_less) = [];
              cand_idx = find(search_series_current==max(search_series_less),1);
              cand_idx_less = find(search_series_less==max(search_series_less),1);
          else
              not_done = false;
          end
      end
      cand_idx_filt = cand_idx-((filt_dist-1)/2):1:cand_idx+((filt_dist-1)/2);
      coeff = prctile(y(cand_idx_filt),pct);
  end
end