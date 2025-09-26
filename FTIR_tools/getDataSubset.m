function [y_subset, x_subset] = getDataSubset(x,y,x_range)
dx = abs(x(2) - x(1));
l = x > x_range(1) - dx/2 & x < x_range(2) + dx/2;
x_subset = x(l);
y_subset = y(l);
end