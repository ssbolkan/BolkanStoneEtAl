function rasterplot(x,y)
h=line(repmat(y,2,1),[x.*ones(1,length(y));(x+.9).*ones(1,length(y))]);
set(h,'Color','k')
hold on
end