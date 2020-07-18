function [X,Y,Z,C] = get_depth_from_viewer(V,F,n,C0)
% Given a surface mesh, a grid size and a current angle from the viewer,
% extract a hexagonal grid simulating what would happen if we did a 3D
% scan of the surface from our current viewpoint.

if nargin==3
    C0 = zeros(size(V,1),3);
end

cam_position = campos;
cam_target = camtarget;
normal = cam_target-cam_position;
normal = normal./norm(normal);

hx = 1/n;
hy = hx*sin(pi/3);
d = hx*cos(pi/3);
[X,Y] = meshgrid(-1.5:hx:1.5,1.5:-hy:-1.5);
rows_X = size(X,1);
assert(mod(rows_X,2)==1);
cols_X = size(X,2);
X(2:2:rows_X-1,:) = X(2:2:rows_X-1,:)-d;
sample_plane = [X(:),Y(:),zeros(size(X(:),1),1)];
if abs(dot(normal',[0;0;1]))<1-1e-5
R = rotmatrix([0;0;1],-normal');
theta = acos(dot(normal',[0;0;1]));
%R = R*[cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1];
else
    R = dot(-normal,[0;0;1]).*eye(3,3);
end
sampling_plane = (R*(sample_plane'))';
sampling_plane = sampling_plane + repmat(cam_position,size(sampling_plane,1),1);

D = repmat(normal,size(sampling_plane,1),1);
[flag, T, lambda] = ray_mesh_intersect(sampling_plane,D, V, F);
Z = reshape(-T,size(X,1),size(X,2));

R = zeros(size(sampling_plane,1),1);
G = zeros(size(sampling_plane,1),1);
B = zeros(size(sampling_plane,1),1);
C = zeros(size(X,1),size(X,2),3);

for i=1:size(sampling_plane,1)
    if flag(i)~=0
        R(i) = lambda(i,1)*C0(F(flag(i),1),1)+lambda(i,2)*C0(F(flag(i),2),1)+...
            lambda(i,3)*C0(F(flag(i),3),1);
        G(i) = lambda(i,1)*C0(F(flag(i),1),2)+lambda(i,2)*C0(F(flag(i),2),2)+...
            lambda(i,3)*C0(F(flag(i),3),2);
        B(i) = lambda(i,1)*C0(F(flag(i),1),3)+lambda(i,2)*C0(F(flag(i),2),3)+...
            lambda(i,3)*C0(F(flag(i),3),3);
    end
end

R = reshape(R,size(X,1),size(X,2));
G = reshape(G,size(X,1),size(X,2));
B = reshape(B,size(X,1),size(X,2));
C(:,:,1) = R;
C(:,:,2) = G;
C(:,:,3) = B;

Z(Z>-Inf) = Z(Z>-Inf)-min(min(Z(Z>-Inf)))+0.05;
Z(Z<0) = 0;



remove_rows = [];
remove_cols = [];

for i=1:size(X,1)
    if i>size(X,1)-5
        if ~any(any(Z((i-3):i,:)))
        remove_rows = [remove_rows,i];
        end
    elseif i<5
    if ~any(any(Z(i:i+4,:)))
        remove_rows = [remove_rows,i];
    end
    else
    if ~any(any(Z((i-4):i+5,:)))
        remove_rows = [remove_rows,i];
    end
    end
end

% for i=1:size(X,2)
%     if i>size(X,2)-5
%         if ~any(any(Z(:,(i-3):i)))
%         remove_cols = [remove_cols,i];
%         end
%     elseif i<5
%     if ~any(any(Z(:,i:i+4)))
%         remove_cols = [remove_cols,i];
%     end
%     else
%     if ~any(any(Z(:,(i-4):i+5)))
%         remove_cols = [remove_cols,i];
%     end
%     end
% end
% 
% 
% Z(remove_rows,:) = [];
% X(remove_rows,:) = [];
% Y(remove_rows,:) = [];
% Z(:,remove_cols) = [];
% X(:,remove_cols) = [];
% Y(:,remove_cols) = [];




end

