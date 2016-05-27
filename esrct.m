%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=esrct(x,t)
% 1MW, 90 fuel elements, diameter = 3.4..cm, fuel height = 38 cm)
% 1e6/90/(pi*.034823^2/4*.381)
y=30.6e6*min(1,t/10);
end
