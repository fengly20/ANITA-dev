function map_mat = valueChangeMat(results_cells, missing_value,low_high)

if nargin<3
    low_high = [2 98];
end

[rows,cols] = size(results_cells);
map_mat_raw = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols
       map_mat_raw(i,j) = results_cells{i,j}{3}(end)-results_cells{i,j}{3}(1) ;
    end
end

map_mat_raw(isnan(map_mat_raw)) = missing_value;

%contrast stretch the ouput image
map_mat = vi_stretch(map_mat_raw,low_high);

end