function dist_mat = calDistance(knot_set,coeff_set,pts)

  for m = 1:(length(knot_set)-1 )
      edge = [knot_set(m) coeff_set(m) knot_set(m+1) coeff_set(m+1)];
      dist_mat(:,m) = distancePointEdge(pts, edge);
  end % end of for m = 1:length(knot_set)-1