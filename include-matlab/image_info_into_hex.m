function [X,Y,Z,hex_edge_density] = image_info_into_hex(IM_INPUT,E)
% Interpolate values from a square grid into a hexagonal grid.

n_input = size(IM_INPUT,2);
m_input = size(IM_INPUT,1);
% hx_input = 1/n_input;
% hy_input = 1/m_input;
[X_input,Y_input] = meshgrid(linspace(-0.1,1.1,n_input),linspace(-0.1,1.1,m_input));
n = size(X_input,2);
hx = 1/n;
hy = hx*sin(pi/3);
d = hx*cos(pi/3);
[X,Y] = meshgrid(-0.1:hx:1.1,1.1:-hy:-0.1);
m = size(X,1);
n = size(X,2);
if mod(m,2)==1
    X(2:2:m-1,:) = X(2:2:m-1,:)-d;
else
    X(2:2:m,:) = X(2:2:m,:)-d;
end
IM_INPUT(IM_INPUT==0) = NaN;
Z = interp2(X_input,Y_input,IM_INPUT,X,Y);
Z(isnan(Z)) = 0;
A = [1,1;1,1];
if nargin>1
E = imdilate(E,A);
hex_edge_density = interp2(X_input,Y_input,E,X,Y);
hex_edge_density(isnan(hex_edge_density)) = 0;
%hex_edge_density = 0;
end
end
