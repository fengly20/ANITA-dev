function [R_sub] = geoRefSubset(R_ori,xmin,xmax,ymin,ymax)

    if nargin<5
        error('need all variables, chump')
    end
    
    old_RasterSize = R_ori.RasterSize;
    new_RasterSize = [ymax-ymin+1 xmax-xmin+1];
    
    old_XWorldLimits = R_ori.XWorldLimits;
    old_YWorldLimits = R_ori.YWorldLimits;
    
    new_XWorldLimits = [old_XWorldLimits(1)+(xmin-1)*30 old_XWorldLimits(2)-(old_RasterSize(2)-xmax)*30];
    new_YWorldLimits = [old_YWorldLimits(1)+(ymin-1)*30 old_YWorldLimits(2)-(old_RasterSize(1)-ymax)*30];
      
    %set properties of spatial reference
    R_sub = R_ori;
    R_sub.XWorldLimits = new_XWorldLimits;
    R_sub.YWorldLimits = new_YWorldLimits;
    R_sub.RasterSize = new_RasterSize;
   
end