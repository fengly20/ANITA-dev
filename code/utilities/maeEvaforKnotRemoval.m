function mae_iter_ortho = maeEvaforKnotRemoval(knot_set,coeff_set,pts,y_pos_idx, pct)

for i=1:length(knot_set)-2 % evaluating only inner knots so the knot index would be i+1 as below
    remove_idx = i+1;
    include_idx = setdiff(1:length(knot_set),remove_idx);
    new_knot_set = knot_set(include_idx);
    new_coeff_set = coeff_set(include_idx);
  
    dist = calDistance(new_knot_set,new_coeff_set,pts);
                            
    ortho_err = min(dist,[],2);
    ortho_err(y_pos_idx)=ortho_err(y_pos_idx)*pct;
    ortho_err(~y_pos_idx) = ortho_err(~y_pos_idx)*(100-pct);
    mae_iter_ortho(i+1) = mean(ortho_err);
end
