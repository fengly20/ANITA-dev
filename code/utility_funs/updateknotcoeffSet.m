function [knot_set,coeff_set,loc_set] = updateknotcoeffSet(knot_set,coeff_set,loc_set,x,cand_loc,coeff)
knot_loc = x(cand_loc);
knot_set = sort([knot_set; knot_loc]);
new_loc_idx = find(knot_set==knot_loc);
coeff_set = [coeff_set(1:new_loc_idx-1);coeff;coeff_set(new_loc_idx:end)];
loc_set = sort([loc_set; cand_loc]);
end

  