clear; close all % start clean

% set up some random stupid geometry, X,Y are control points
X = [ 0, 1, 2, 3, 3;
     -1, 0, 1, 2, 2;
      0, 1, 2, 2, 3;
     -1, 0, 1, 1, 1;
     -2,-1,-1,-1,-1;
     -3,-2,-2,-2,-2];
Y = [0, 0, 0, 0, 1;
     1, 1, 1, 1, 2;
     2, 2, 3, 4, 5;
     3, 3, 4, 5, 6;
     4, 5, 6, 7, 8;
     4, 5, 6, 7, 8];
knot1 = [0,0,0,1,2,3,3,3];
knot2 = [0,0,0,1,2,3,4,4,4];


% compute basis function for all evaluation points
xi = linspace(knot1(1), knot1(end), 50);
eta= linspace(knot2(1), knot2(end), 60);
[Nu, dNu] = getBSplineBasisAndDerivative(2, xi,  knot1);    % evaluation mesh
[Nku, ~]  = getBSplineBasisAndDerivative(2, knot1,  knot1); % just used for higlighting element boundaries
ddNu      = getBSplineHighDerivative(    2, xi,  knot1, 2); % double derivative for piola mapping

[Nv, dNv] = getBSplineBasisAndDerivative(2, eta, knot2);    % and all the same just for eta-values
[Nkv, ~]  = getBSplineBasisAndDerivative(2, knot2,  knot2);
ddNv      = getBSplineHighDerivative(    2, eta, knot2, 2);


% plot the mesh we're looking at
x  = Nu'  * X' * Nv; % evaluation points
y  = Nu'  * Y' * Nv;
xu = Nku' * X' * Nv; % knot lines (element boundaries, u-direction)
yu = Nku' * Y' * Nv;
xv = Nu'  * X' * Nkv; % knot lines (element boundaries, v-direction)
yv = Nu'  * Y' * Nkv;
figure; hold on;
  plot(x, y,  'k-', 'Color', [.6, .6, .6]); % mesh points, light gray
  plot(x',y', 'k-', 'Color', [.6, .6, .6]);
  plot(xv, yv,  'k-', 'LineWidth', 2);      % element lines, fat black
  plot(xu',yu', 'k-', 'LineWidth', 2);
  plot(X, Y,  'bo-');                       % control points, blue dots
  plot(X',Y', 'b-');
  title('Physical mesh $$\mathbf{x}(\xi, \eta)$$', 'Interpreter', 'LaTeX')


% For now we are going to look at the following function (a single basisfunction):
U = [0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0];
V = [0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 0 0;
     0 0 0 1 0;
     0 0 0 0 0];


% Display function using REGULAR mapping (no Piola)
u  = Nu' * U' * Nv; % evaluating vector on the specified points
v  = Nu' * V' * Nv;
dudx = zeros(size(u));
dudy = zeros(size(u));
dvdx = zeros(size(v));
dvdy = zeros(size(v));
figure;
  sgtitle('REGULAR mapping')
  subplot(1,3,1); hold on;
    plot(xv, yv,  'k-', 'LineWidth', 2);      % element lines, fat black
    plot(xu',yu', 'k-', 'LineWidth', 2);
    quiver(x,y,u,v);
    title('Vector basisfunction [u,v]');
  subplot(1,3,2);
    surf(x,y,zeros(size(x)),u); hold on;
    view(2);
    colorbar;
    plot(xv, yv,  'k-', 'LineWidth', 2);      % element lines, fat black
    plot(xu',yu', 'k-', 'LineWidth', 2);
    title('First component u');
  subplot(1,3,3);
    surf(x,y,zeros(size(x)),v); hold on;
    view(2);
    colorbar;
    plot(xv, yv,  'k-', 'LineWidth', 2);      % element lines, fat black
    plot(xu',yu', 'k-', 'LineWidth', 2);
    title('Second component v');

% Display function using PIOLA mapping
dxdxi  = dNu' * X' *  Nv; % compute jacobian everywhere
dxdeta =  Nu' * X' * dNv;
dydxi  = dNu' * Y' *  Nv;
dydeta =  Nu' * Y' * dNv;
cp = [X(:)';Y(:)'];
for i=1:numel(xi)
  for j=1:numel(eta)
    N = [kron(  Nu(:,i),  Nv(:,j))';
         kron( dNu(:,i),  Nv(:,j))';
         kron(  Nu(:,i), dNv(:,j))';
         kron( dNu(:,i), dNv(:,j))';
         kron(ddNu(:,i),  Nv(:,j))';
         kron( dNu(:,i), dNv(:,j))';
         kron(  Nu(:,i),ddNv(:,j))'];
    map = computeGeometry(cp,N);

    u_parametric = [Nu(:,i)' * U' * Nv(:,j);
                    Nu(:,i)' * V' * Nv(:,j)];
    u_parametric_diff = [dNu(:,i)' * U' * Nv(:,j);
                            dNu(:,i)' * V' * Nv(:,j);
                            Nu(:,i)' * U' * dNv(:,j);
                            Nu(:,i)' * V' * dNv(:,j)];
    u_physical = 1/map.detJ * map.J * u_parametric;
    [u_physical,u_physical_derivative] = piolaTransform(map,u_parametric,u_parametric_diff);
    u(i,j) = u_physical(1);
    v(i,j) = u_physical(2);
    dudx(i,j) = u_physical_derivative(1);
    dvdx(i,j) = u_physical_derivative(2);
    dudy(i,j) = u_physical_derivative(3);
    dvdy(i,j) = u_physical_derivative(4);
  end
end
figure;
  sgtitle('PIOLA mapping')
%   subplot(3,3,1); hold on;
%     plot(xv, yv,  'k-', 'LineWidth', 2);      % element lines, fat black
%     plot(xu',yu', 'k-', 'LineWidth', 2);
%     quiver(x,y,u,v);
%     title('Vector basisfunction [u,v]');
G = {u,dudx,dudy,v,dvdx,dvdy};
T = {'First component u','x-derivative of first component u','y-derivative of first component u',...
    'Second component v','x-derivative of second component v','y-derivative of second component v'};
Gmax = zeros(1,6);
Gmin = zeros(1,6);
for i = 1:6
    Gmax(i) = max(max(G{i}));
    Gmin(i) = min(min(G{i}));
end
Gmax = max(Gmax);
Gmin = min(Gmin);
w = winter;
a = autumn;
cm = [flip(w); a];
L = max(abs(Gmax),abs(Gmin));
for i = 1:6
    subplot(2,3,i);
    surf(x,y,zeros(size(x)),G{i}); hold on;
    view(2);
    colorbar;
    if ~((i == 1)||(i == 4))
        clim([-1.4,1.4]);
        colormap(cm)
    end
    plot(xv, yv,  'k-', 'LineWidth', 2);
    plot(xu',yu', 'k-', 'LineWidth', 2);
    title(T{i});
end
%   subplot(2,3,1);
%     surf(x,y,zeros(size(x)),u); hold on;
%     view(2);
%     colorbar;
%     plot(xv, yv,  'k-', 'LineWidth', 2);
%     plot(xu',yu', 'k-', 'LineWidth', 2);
%     title('First component u');
%   subplot(2,3,2);
%     surf(x,y,zeros(size(x)),dudx); hold on;
%     view(2);
%     colorbar;
%     plot(xv, yv,  'k-', 'LineWidth', 2);
%     plot(xu',yu', 'k-', 'LineWidth', 2);
%     title('x-derivative of first component u');
%   subplot(2,3,3);
%     surf(x,y,zeros(size(x)),dudy); hold on;
%     view(2);
%     colorbar;
%     plot(xv, yv,  'k-', 'LineWidth', 2);
%     plot(xu',yu', 'k-', 'LineWidth', 2);
%     title('y-derivative of first component u');
%   subplot(2,3,4);
%     surf(x,y,zeros(size(x)),v); hold on;
%     view(2);
%     colorbar;
%     plot(xv, yv,  'k-', 'LineWidth', 2);
%     plot(xu',yu', 'k-', 'LineWidth', 2);
%     title('Second component v');
%   subplot(2,3,5);
%     surf(x,y,zeros(size(x)),dvdx); hold on;
%     view(2);
%     colorbar;
%     plot(xv, yv,  'k-', 'LineWidth', 2);
%     plot(xu',yu', 'k-', 'LineWidth', 2);
%     title('x-derivative of second component v');
%   subplot(2,3,6);
%     surf(x,y,zeros(size(x)),dvdy); hold on;
%     view(2);
%     colorbar;
%     plot(xv, yv,  'k-', 'LineWidth', 2);
%     plot(xu',yu', 'k-', 'LineWidth', 2);
%     title('y-derivative of second component v');
