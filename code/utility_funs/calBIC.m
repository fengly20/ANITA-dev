function bic_remove = calBIC(ortho_err,keep_knots,penalty)
% BIC acccumulation
 
  positive_idx = ortho_err>0; %in case a value is exactly 0
% this line is maybe the most computationally expensive single line in nita
  dist_ortho_err = fitdist(ortho_err(positive_idx), 'lognormal');
  loglik = -dist_ortho_err.negloglik; % log likelyhood
  num_segs=length(keep_knots)-1;
%user-defined penalty parameter
  bic_remove = -2*loglik+penalty*num_segs*log(length(ortho_err));
end
