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
  
  if sum(mov_std) == 0
      mov_std = mov_std+1;
  end
  
  mov_cv = mov_mean./mov_std;        
  
  if nargin > 5
      zscore = (x-mov_mean)./mov_std;
      search_series = abs(zscore);
  else 
      search_series = mov_cv;
  end

  % in case that search_series has two (or more) max values equal to each
  % other, add a bit white noise to ensure that each value is unique 
  if length(unique(search_series))~=length(search_series)
      noise_amplitude = 0.00000000000001;
      noise = noise_amplitude*randn(1, length(search_series))';
      search_series = search_series+noise;
  end
  
  % candidate point needs to be found out of filt_dist of beginning or end
  % of time series
  search_series_inner = search_series((filt_dist+1):(length(search_series)-filt_dist));
  cand_idx = find(search_series==max(search_series_inner),1);
  
  cand_idx_filt = cand_idx-((filt_dist-1)/2):1:cand_idx+((filt_dist-1)/2);
  coeff = prctile(y(cand_idx_filt),pct);
end 
