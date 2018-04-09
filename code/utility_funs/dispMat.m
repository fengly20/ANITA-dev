function dummy = dispMat(mat)

stretched_mat = vi_stretch(mat,[1 99]);
figure, img = imagesc(stretched_mat); axis equal; colorbar
set(img,'AlphaData',~isnan(mat))
end

