function [knot_set,coeff_set,coeff_indices] = updateknotcoeffSet(knot_set,coeff_set,coeff_indices,knot_loc,cand_idx,coeff)
  knot_set = sort([knot_set; knot_loc]);
  new_loc_idx = find(knot_set==knot_loc);
  coeff_set = [coeff_set(1:new_loc_idx-1);coeff;coeff_set(new_loc_idx:end)];
  coeff_indices = sort([coeff_indices; cand_idx]);
  