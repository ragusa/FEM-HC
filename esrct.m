function y=esrct(x,t)
% 1MW, 90 fuel elements, Zr inner rod radius=0.003175m
% fuel meat radius=0.0174115m, fuel height = 0.381m)
% 1e6/90/(pi*(0.0174115^2-0.003175^2)*.381)
y=3.1674e7*min(1,t/10);
end