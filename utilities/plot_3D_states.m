function [x,y] = plot_3D_states(filename,emin,emax,num,broadmeV)
% ---------------------------------------------
% matlab function to plot 3D states broadened by
% a lorentizan
% ---------------------------------------------
X = load(filename);

n = num;
x = linspace(emin,emax, n);
y = zeros(n,1);

g = broadmeV;

for ii = 1:n
    for vv = 1:length(X)
        y(ii) = y(ii) + (pi*g/2)*(1.0/pi)*(g/2)/((x(ii) - X(vv))^2 + (g/2)^2);
    end
end



