function [img_yeardate,img_doy] = landsatImgDate(img_name) 
L_idx = strfind(img_name, 'L');
img_yeardate = img_name(L_idx+12:L_idx+19);
yymmdd = [str2num(img_yeardate(1:4)), str2num(img_yeardate(5:6)), str2num(img_yeardate(7:8))]; 
img_yeardate = str2num(img_yeardate);
img_doy = day(datetime(yymmdd), 'dayofyear');
end


