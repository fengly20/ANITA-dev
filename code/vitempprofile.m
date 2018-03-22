function [noise, include_idx, y]=vitempprofile(veg_index,year_vec,knot_set,coeff_set,xi,yi,axis_in, c)

%this function is only for plotting.
%inputs
%veg_index: the x by y by num_dates index image
%year_vec: vector of only the years distilled from date_num
%knot_set: the cell array of all knots processed from the veg_index
%coeff_set: same as knot_set but for coefficients
%xi, yi: the pixel of interest
%axis_in: 1x4 vector of axis extents [min(x) max(x) min(y) max(y)]
%vecort of colors for dots (e.g. image dy of year

%check if input is an image or a single pixel
if size(veg_index,2)>1
    y = squeeze(veg_index(yi,xi,:));
else
    y = veg_index;
end

%put in an index that will change depending on dataset
%include_idx = find(y<100 & not(y==0) & y>-100);
%include_idx = find(y>-100 & y<100);
%include_idx = find(y>0 & y<10000);
include_idx = find(not(isnan(y)));


%if user provides a "c" for scatter coloring do the scatter plot (bottom)
if nargin<8
    figure, plot(year_vec(include_idx),y(include_idx),'.'); hold on, box on
        knots_x = knot_set{yi}{xi};
        %knots_x_fulldate = datetime(knot_set{yi}{xi},'ConvertFrom','juliandate');
        %knots_x = datevec(knots_x_fulldate);
        coeffs_y = coeff_set{yi}{xi};
        plot(knots_x(:,1), coeffs_y,'ro')
        plot(knots_x(:,1), coeffs_y,'-r')
        %optionally allow user to set axes
        if nargin==7
            axis(axis_in)
        end
        hold off
else
    figure, scatter(year_vec(include_idx),y(include_idx), 20, c(include_idx), 'filled'); hold on, box on
        
        if iscell(knot_set)
            knots_x = knot_set{yi}{xi};
            coeffs_y = coeff_set{yi}{xi};
        else
            knots_x = knot_set;
            coeffs_y = coeff_set;
        end
        plot(knots_x(:,1), coeffs_y,'ro')
        plot(knots_x(:,1), coeffs_y,'-r')
        axis(axis_in)
        colorbar;
        
        hold off    
end

noise = mean(abs(diff(y(include_idx))));


end
