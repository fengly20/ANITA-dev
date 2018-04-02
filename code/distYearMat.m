function map_mat = distYearMat(metrics_cells,missing_value,low_high) 

if nargin<3
    low_high = [2 98];
end

[rows,cols] = size(metrics_cells);
map_mat_raw = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols
       map_mat_raw(i,j) = (metrics_cells{i,j}{3}+metrics_cells{i,j}{4})/2;
    end
end

map_mat_raw(isnan(map_mat_raw)) = missing_value;

map_mat = vi_stretch(map_mat_raw,low_high);
end