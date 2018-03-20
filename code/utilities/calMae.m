function mae = calMae(dist)
dist = min(dist,[],2);
mae = mean(dist); % mean or median   
end