clf
n = 1000;
hx = 1/n;
hy = hx*sin(pi/3);
d = hx*cos(pi/3);
[X,Y] = meshgrid(-0.1:hx:1.1,1.1:-hy:-0.1);
m = size(X,1);
n = size(X,2);
X(2:2:m-1,:) = X(2:2:m-1,:)-d;
Z0 = .5.*ones(size(X,1),size(X,2));

% PART THAT IS ZERO
b = find((Y(:)<0.2+(X(:)-.5).^2 ));
bc = zeros(size(b,1),1);
b = find((Y(:)<0.2+0.2.*cos(10.*X(:))));
bc = zeros(size(b,1),1);



% PART THAT IS ONE
b2 = find((Y(:)>0.6+(X(:)-.5).^2 ));
bc2 = ones(size(b2,1),1);
b2 = find((Y(:)>0.6+0.2.*sin(10.*X(:))));
bc2 = ones(size(b2,1),1);


% Together
b = [b;b2];
bc = [bc;bc2];
[b,IA,IC] = unique(b);
bc = bc(IA);


Z0(b) = bc;
Z0(find((X(:)>1)+(X(:)<0)+(Y(:)>1)+(Y(:)<0))) = 0;

Z = sparsify_height_field_admm(X,Y,Z0,'InterpolateB',b,'InterpolateBC',bc);
 %Z = sparsify_height_field_admm(X,Y,Z0);
 
 
 
 clf
 system('mkdir ribbon'); save_everything(X,Y,Z,Z0,'ribbon',false);
 
 [V,F] = read_triangle_mesh('ribbon/objs/output.obj');
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
 
 Z0_plot = Z0;
 Z0_plot(setdiff(1:size(Z0,1)*size(Z0,2),b)) = NaN;
 %Z0_plot(find(((X(:)-.5).^2 + (Y(:)-.5).^2)>(0.305.^2))) = NaN;
 surf(X-1.2,Y,Z0_plot,'EdgeColor','none')
 
 colormap(cbrewer('RdYlBu',20))

 axis equal
 axis off
 grid off
 view([0 26])
 camlight
 set(gcf,'Color','w');
  add_isolines()
 drawnow
 %imwrite(myaa({'raw',5}),'ribbon.png')
 
 
 
 
 