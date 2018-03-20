function map_mat = dateThresMat(results_cells,thres,missing_value)

[rows,cols] = size(results_cells);
map_mat = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols
        first_knot = results_cells{i,j}{2}(1);
        last_knot = results_cells{i,j}{2}(end);
        first_coeff = results_cells{i,j}{3}(1);
        last_coeff = results_cells{i,j}{3}(end);
        
        x = [first_knot last_knot];
        y = [first_coeff last_coeff];
        
        ab_set = [[1; 1] x(:)]\y(:); % Calculate Parameter Vector
        slope = ab_set(2);
        intercept = ab_set(1);
        
        potential_knot = (thres-intercept)/slope; 
        
        if potential_knot < first_knot
            potential_knot = NaN;
        elseif potential_knot > last_knot
            potential_knot = NaN;
        else 
            potential_knot = potential_knot; 
        end
            
        map_mat(i,j) = potential_knot;
    end
end
map_mat(isnan(map_mat)) = missing_value;
end