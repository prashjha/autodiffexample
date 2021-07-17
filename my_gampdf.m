function [f] = my_gampdf(z, a, b)

% my_gampdf
% Implementation of gampdf that is compatible with auto differentiation
% in Matlab.
% 
% Takes same inputs as Matlab's gampdf.
%
% Author: Adrian Celaya (github: aecelaya)
% Last modified: 06.28.2021

% Direct implementation of gamma pdf. Does not allow parameters a and b to
% be optimization variables. Optimization toolbox does not like b .^ a.
%f = (1 ./ ((b .^ a) .* my_gamma(a))) .* z .^ (a - 1) .* exp(-z ./ b);

% Use exp(log(gampdf)) to get rid of b .^ a. 
f = exp(((a - 1) .* log(z) - (z ./ b)) - (a .* log(b) + log(my_gamma(a))));

return 