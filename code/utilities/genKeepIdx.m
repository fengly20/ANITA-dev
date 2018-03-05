function keep_idx = genKeepIdx(keep_knots,keep_coeffs,pts,pct,y_pos_idx)

if length(keep_knots)-2 > 0         
  %isolate data points above interp fit.                      
    mae_iter_ortho = maeEvaforKnotRemoval(keep_knots,keep_coeffs,pts,y_pos_idx,pct);
  %now find index for actual knot removal based on
  %minimizing the mae_iter_ortho error
    remove_idx = find(mae_iter_ortho == min(mae_iter_ortho(2:end)),1);
  %keep all knots EXCEPT for the remove_idx
    keep_idx = setdiff(1:length(keep_knots),remove_idx);
else
  %this will only be first and last knots
    keep_idx = [1 3];
end % end of length(keep_knots)-2 > 0