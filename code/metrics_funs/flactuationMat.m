function map_mat = flactuationMat(date,results_cells,metrics_cells,cum_rise_thres)

value_mat = valueChangeMat(-9999,date,'value',metrics_cells);
normrise_mat = normRiseMat(results_cells);

map_mat = value_mat;
map_mat(normrise_mat<cum_rise_thres) = NaN;

end


