function [vi_adj] = vi_stretch(vi,low_high)

%inputs
%Any image. The purpose of this function is to work beyond just intensity
%images (stretchlim).
%[low high] is entered as percentiles. Default is [2 98] indicating a 2%
%clip.

if nargin<2
    low_high = [2 98];
end

low_per = prctile(vi(:), low_high(1));
high_per = prctile(vi(:), low_high(2));

%out of bounds binaries
high_bin = vi>high_per;
low_bin = vi<low_per; 

%swap in saturate value
vi_adj = vi;
vi_adj(high_bin) = high_per;
vi_adj(low_bin) = low_per;




end