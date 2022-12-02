function [N dN] = piolaTransform(map, N_in, dN_in)

N = map.J * N_in / map.detJ;
dN = zeros(4,size(N,2));




% REPLACE THIS ABDULLAH!
% the output matrix dN should contain 4 rows in order: u_1,1  u_2,1  u_1,2  u_2,2 (einstein notation)
% and it should have as many columns as there are basis functions, i.e. check the size of N
