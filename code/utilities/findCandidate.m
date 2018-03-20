% findCandidate
% input: 
% return: non_nan_idx - a n*1 list taht contains indeces for good data
% points.
% note: a bit discussion about the sign for dist, mov_cv, and zscore:
%       The distances are calculated as orthogonal distance between points
%       to an edge, so the dists should only contain positive number. Thus
%       the mov_cv(see below in function for calculation) should also
%       contain only positive numbers. On the other hand, 
%       the zscore(see below in function for calculation) may contain
%       posistive numbers and negative numbers. For breakpoint detection,
%       the large chnage point can be identified as a breakpoint no matter
%       it's posistive or negative so abs() is applied on zscore. 
function [cand_idx,coeff,search_series] = findCandidate(dist,filt_dist,pct,y,x)

  dist = min(dist,[],2);

  mov_mean = movmean(dist,filt_dist);
  mov_std = movstd(dist,filt_dist);
  
  mov_cv = mov_mean./mov_std;        
  
  if nargin > 5
      zscore = (x-mov_mean)./mov_std;
      search_series = abs(zscore);
  else 
      search_series = mov_cv;
  end
      
  cand_idx = find(search_series==max(search_series),1);
    
% if the candidate point is within filt_dist of beginning or end
% of time series, ignore and search for another point. 
% there's a possiblity of the next found cand_idx is also with the filt_dist
% of beginning or end of time series so use while loop to ensure 

  not_done = true;
  search_series_less = search_series;
  cand_idx_in_ssless = cand_idx;
  while not_done
      if cand_idx>(length(y)-filt_dist) || cand_idx<(filt_dist+1)
          search_series_less(cand_idx_in_ssless) = [];
          cand_idx = find(search_series==max(search_series_less),1);
          cand_idx_in_ssless = find(search_series_less==search_series(cand_idx));
      else
          not_done = false;
      end
  end
  
  cand_idx_filt = cand_idx-((filt_dist-1)/2):1:cand_idx+((filt_dist-1)/2);
  coeff = prctile(y(cand_idx_filt),pct);
end 
