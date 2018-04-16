function norm_rise = normRiseMat(results_cells) 

[rows,cols] = size(results_cells);

norm_rise = zeros(rows,cols);
for i = 1:rows
    for j = 1:cols
        results_cell = results_cells{i,j};
        sum_of_abs_rises = sum(abs(results_cell{7}));
        complexity = results_cell{1};
        %norm_rise(i,j) = sum_of_abs_rises/complexity;
        norm_rise(i,j) = sum_of_abs_rises;
    end
end

end