function map = computeGeometry(cp, N)

% The main quantities of the geometry mapping.
x = cp * N(1,:)'  ;                     % Physical coordinate point (x,y)
J = cp * N(2:3,:)';                     % Jacobian matrix [dx/du,dx/dv; dy/du, dy/dv]
H = cp * N(4:6,:)';                     % Hessian matrix [dx/du2, dx/dudv, dx/dv2; dy/du2, dy/dudv, dy/dv2]
H = [H(:,1), H(:,2), H(:,2), H(:,3)];
H = reshape(H,2,2,2);                   % Hessian matrix H(i,j,k)=d^2x(i)/dxi(j)/dxi(k)
invJ = inv(J);                          % Inverse Jacobian [du/dx, du/dy; dv/dx, dv/dy]
detJ = det(J);                          % Jacobi-determinant

% First-order derivative facilities.
Jt1 = H(1,1,1)*J(2,2)+J(1,1)*H(2,1,2)-H(1,1,2)*J(2,1)-J(1,2)*H(2,1,1);
Jt2 = H(1,1,2)*J(2,2)+J(1,1)*H(2,2,2)-H(1,2,2)*J(2,1)-J(1,2)*H(2,1,2);
Jx = invJ(1,1)*[H(1,1,1),H(1,1,2);H(2,1,1),H(2,1,2)]+invJ(2,1)*[H(1,1,2),H(1,2,2);H(2,1,2),H(2,2,2)];
Jy = invJ(1,2)*[H(1,1,1),H(1,1,2);H(2,1,1),H(2,1,2)]+invJ(2,2)*[H(1,1,2),H(1,2,2);H(2,1,2),H(2,2,2)];
detJx = invJ(1,1)*Jt1+invJ(2,1)*Jt2;
detJy = invJ(1,2)*Jt1+invJ(2,2)*Jt2;

map = struct('x',x,'J',J,'H',H,'detJ',detJ,'invJ',invJ,'Jx',Jx,'Jy',Jy,'detJx',detJx,'detJy',detJy);
