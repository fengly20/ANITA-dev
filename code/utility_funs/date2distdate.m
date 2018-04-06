function dist_date = date2distdate(YYYYMMDD)
date_str = num2str(YYYYMMDD);
date_str_format = [date_str(1:4) '-' date_str(5:6) '-' date_str(7:8)];
d = datetime(date_str_format);
doy = day(d,'dayofyear');
dist_date = str2num(date_str(1:4))*1000+doy/365*1000;
end