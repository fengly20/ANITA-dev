function date_num = genDisDate(img_yeardate,img_doy)

if size(img_yeardate,2)>size(img_yeardate,1)
    img_yeardate = img_yeardate';
end
if size(img_doy,2)>size(img_doy,1)
    img_doy = img_doy';
end

date_str = num2str(img_yeardate);
year_vec = str2num(date_str(:,1:4));
partial_year = (img_doy/365)*1000;
date_num = (year_vec*1000)+round(partial_year);

end
