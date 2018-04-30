function map_mat = noiseMat(results_cells)

[rows,cols] = size(results_cells);
map_mat = zeros(rows,cols);

for i = 1:rows
    for j = 1:cols
       map_mat(i,j) = results_cells{i,j}{6};
    end
end

map_mat(map_mat==-999) = NaN;

end
