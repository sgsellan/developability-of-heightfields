% clear;
% [V,F] = readOBJ('church.obj');
% [X, Y, Z] = mesh_to_height_field(V, F, 800, true);
% surf(X,Y,Z)
% surf(X,Y,Z,'EdgeColor','none')
% RX = @(theta) [1,0,0;0,cos(theta), sin(theta);0,-sin(theta),cos(theta)];
% V = V*RX(pi);
% V = V*RX(pi/2);
% [X, Y, Z] = mesh_to_height_field(V, F, 800, true);
% surf(X,Y,Z,'EdgeColor','none')
% RY = @(theta) [cos(theta), 0, sin(theta);0,1,0;-sin(theta),0,cos(theta)];
% RZ = @(theta) [cos(theta), sin(theta),0;-sin(theta),cos(theta),0;...
% 0,0,1];
% V = V*RZ(pi/2);
% [X, Y, Z] = mesh_to_height_field(V, F, 1000, true);
% surf(X,Y,Z,'EdgeColor','none')
% save('church_above.mat','X','Y','Z');

A = imread('test.png');
A = A(1:1:end,1:1:end,1);
if size(A,1)>size(A,2)
    A = A((size(A,1)-size(A,2)+1):size(A,1),1:size(A,2));
else
    A = A(1:size(A,1),1:size(A,1));
end
%A = 255-A;
A(A<250) = 0;
[X,Y,Z0] = image_info_into_hex(double(A));
%Z0 = 255 - Z0;
b = find(Z0(:)>1);

Z0 = Z0./255;

bc = zeros(size(b,1),1);



b2 = find((Y(:)<0.67).*(Y(:)>0.66).*(X(:)>0.1).*(X(:)<0.85));
bc2 = ones(size(b2,1),1);
b3 = find((Y(:)<0.83).*(Y(:)>0.46).*(X(:)>0.38).*(X(:)<0.40));
bc3 = ones(size(b3,1),1);
b4 = find((Y(:)<0.91).*(Y(:)>0.665).*(X(:)>0.78).*(X(:)<0.79));
bc4 = ones(size(b4,1),1);
b2 = unique([b2;b3;b4]);
bc2 = .3.*(-2.*abs(Y(b2)-0.665)+1);

% 
% 
% load('church_above.mat')
% fixed_zero = find(Z(:)>1.3);
% Z0 = ones(size(X,1),size(X,2));
% b = fixed_zero;
% bc = zeros(size(b,1),1);

% 
b = [b;b2];
bc = [bc;bc2];

Z0 = .3.*ones(size(X,1),size(X,2));
Z0(b) = bc;
Z0(find((X(:)>1.07)+(X(:)<-0.07)+(Y(:)>1)+(Y(:)<0.3))) = 0;

Z = sparsify_height_field_admm(X,Y,Z0,'InterpolateB',b,'InterpolateBC',bc);

 
 
 clf
 % Z(find(((X(:)-.5).^2 + (Y(:)-.5).^2)>(0.305.^2))) = NaN;

 system('mkdir church'); save_everything(X,Y,Z,Z0,'church',false);
 
 [V,F] = read_triangle_mesh('church/objs/output.obj');
 hold off
 tsurf(F,V,fsoft,fphong,falpha(1,0));
 O = unique(outline(F));
 theta = linspace(0,2*pi,100);
 theta = theta(1:end);
 U = [cos(theta'),sin(theta'),zeros(size(theta,2),1)];
 bigU = [cos(theta'),sin(theta'),zeros(size(theta,2),1)];
 bigU = bigU.*0.3;
 bigU = bigU+[0.5,0.5,0];
 U = U.*0.1;
 U = U+[0.5,0.5,0.5];
 G = delaunay(U(:,1),U(:,2));
 hold on
 %tsurf(G,U-[0.8,0,0],fsoft,fphong,falpha(1,0));
 %plot3(bigU(:,1)-0.8,bigU(:,2),zeros(size(U,1),1),'Color',[0.6471,0,0.1490])
 
%  Z0_plot = Z0;
%  Z0_plot(setdiff(1:size(Z0,1)*size(Z0,2),b)) = NaN;
%  Z0_plot(find(((X(:)-.5).^2 + (Y(:)-.5).^2)>(0.305.^2))) = NaN;
%  surf(X-1.2,Y,Z0_plot,'EdgeColor','none')
 
 colormap(cbrewer('RdYlBu',20))

 axis equal
 axis off
 grid off
 view([0 26])
 camlight
 set(gcf,'Color','w');
 % add_isolines()
 drawnow
 %imwrite(myaa({'raw',5}),'church.png')
 
 
 
 
 