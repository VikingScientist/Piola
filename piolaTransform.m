function [N,dN] = piolaTransform(map,N_in,dN_in)

% The Piola-mapped function.
N = (map.J*N_in)/map.detJ;

% The Piola-mapped first-order derivatives.
P0 = map.J/map.detJ;
Px = -(map.detJx/map.detJ^2)*map.J+map.Jx/map.detJ;
Py = -(map.detJy/map.detJ^2)*map.J+map.Jy/map.detJ;
dN = zeros(4,size(N,2));
r1 = dN_in(1,:)*map.invJ(1,1)+dN_in(3,:)*map.invJ(2,1);
r2 = dN_in(2,:)*map.invJ(1,1)+dN_in(4,:)*map.invJ(2,1);
r3 = dN_in(1,:)*map.invJ(1,2)+dN_in(3,:)*map.invJ(2,2);
r4 = dN_in(2,:)*map.invJ(1,2)+dN_in(4,:)*map.invJ(2,2);
dN(1:2,:) = Px*N_in+P0*[r1;r2];
dN(3:4,:) = Py*N_in+P0*[r3;r4];