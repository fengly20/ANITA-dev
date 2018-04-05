function map_mat = errorNoiseRatioMat(results_cells)

[rows,cols] = size(results_cells);
map_mat = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols
        linear_eeror_mat(i,j) = results_cells{i,j}{4};
        noise_mat(i,j) = results_cells{i,j}{6};
    end
end

linear_eeror_mat(linear_eeror_mat==-999) = NaN;
noise_mat(noise_mat==-999) = NaN;

map_mat = linear_eeror_mat./noise_mat;

end
