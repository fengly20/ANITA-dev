function [cand_idx,coeff] = findNextCandidate(coeff_indices,search_series_current,filt_dist,pct,y)
  
  search_series_new = search_series_current;
  % remove current selected candidate locations using coeff_indices
  search_series_new(coeff_indices) = [];
  
  % now examine the search_series_new. 
  % if search_series_new == (filt_dist-1)*2, meaning rest locations are all
  % within filt_dist (both begining and the end), so -999 should be
  % assigned to cand_idx and coeff;
  % else select next potential candicate 
  if length(search_series_new) == (filt_dist-1)*2 
      cand_idx = -999;
      coeff = -999;
  else
      search_series_new_inner = search_series_new((filt_dist-1+1):(length(search_series_new)-(filt_dist-1)));
      cand_idx = find(search_series_current==max(search_series_new_inner),1);
      
      cand_idx_filt = cand_idx-((filt_dist-1)/2):1:cand_idx+((filt_dist-1)/2);
      coeff = prctile(y(cand_idx_filt),pct);
  end
end