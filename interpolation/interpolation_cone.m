% n = 1000;
% hx = 1/n;
% hy = hx*sin(pi/3);
% d = hx*cos(pi/3);
% [X,Y] = meshgrid(-0.1:hx:1.1,1.1:-hy:-0.1);
% m = size(X,1);
% n = size(X,2);
% X(2:2:m-1,:) = X(2:2:m-1,:)-d;
% Z0 = ones(size(X,1),size(X,2));
% 
% b = find(((X(:)-.5).^2 + (Y(:)-.5).^2)>(0.3.^2));
% bc = zeros(size(b,1),1);
% b2 = find(((X(:)-.5).^2 + (Y(:)-.5).^2)<(0.1.^2));
% bc2 = ones(size(b2,1),1);
% b = [b;b2];
% bc = [bc;bc2];
% Z0(b) = bc;
% 
% Z = sparsify_height_field_admm(X,Y,Z0,'InterpolateB',b,'InterpolateBC',bc);
%  %Z = sparsify_height_field_admm(X,Y,Z0);
%  
%  
%  
%  
%  system('mkdir cone'); save_everything(X,Y,Z,Z0,'cone');
%  
 [V,F] = read_triangle_mesh('cone/objs/output.obj');
 hold off
 tsurf(F,V.*[1,1,0.5],fsoft,fphong,falpha(1,0));
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
 t = tsurf(G,U-[0.8,0,0],fsoft,fphong,falpha(1,0));
 plot3(bigU(:,1)-0.8,bigU(:,2),zeros(size(U,1),1),'Color',[0.6471,0,0.1490])
 colormap(cbrewer('RdYlBu',20))

 axis equal
 axis off
 grid off
 view([0 26])
 camlight
 set(gcf,'Color','w');
  add_isolines()
 drawnow
 imwrite(myaa({'raw',5}),'cone.png')
 
 
 
 
 