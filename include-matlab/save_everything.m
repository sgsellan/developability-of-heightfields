function [J,J0] = save_everything(X,Y,Z,Z0,directory,boundary)
% Given output and input from sparsify_heightfield_admm, we save input and output
% as triangulated .objs and saves data about principal curvatures, Gaussian curvature,
% etc. and many other quantities needed for debugging and experimental evaluation.

if nargin==5
    boundary = true;
end

system(['mkdir ',directory,'/objs']);
system(['mkdir ',directory,'/pngs']);
save([directory,'/just_depths.mat'],'X','Y','Z','Z0')
m = size(X,1);
n = size(X,2);
always_II = [];
for j=3:(n-2) % for sub2ind but let's do one thing at a time.
    always_II = [always_II, ((j-1)*m)+(3:(m-2))];
end
background = unique([find(Z0(:)==0);find(isnan(Z0(:)));find(Z0(:)==-Inf)]);
II = setdiff(1:m*n,background);
II = intersect(II,always_II);
[NN,BXX,BXY,BYX,BYY] = get_connectivity(X,Y,Z,II,'hex');
V = [X(:),Y(:),Z(:)];
V0 = [X(:),Y(:),Z0(:)];
distance = (Z(:)-Z0(:)).^2;



F = [NN(2,:)',NN(1,:)',NN(4,:)';...
    NN(5,:)',NN(2,:)',NN(4,:)';...
    NN(7,:)',NN(5,:)',NN(4,:)';...
    NN(6,:)',NN(7,:)',NN(4,:)';...
    NN(3,:)',NN(6,:)',NN(4,:)';...
    NN(1,:)',NN(3,:)',NN(4,:)'];
F = remove_duplicate_simplices(F);
F0 = F;
are_non_zero = (V(F(:,1),3)>0).*(V(F(:,2),3)>0).*(V(F(:,3),3)>0);
F(~are_non_zero,:) = [];
are_non_zero = (V0(F0(:,1),3)>0).*(V0(F0(:,2),3)>0).*(V0(F0(:,3),3)>0);
F0(~are_non_zero,:) = [];
[V,I,J] = remove_unreferenced(V,F);
F = I(F);
[V0,I0,J0] = remove_unreferenced(V0,F0);
F0 = I0(F0);
off = [1.4*(max(V(F,1))-min(V(F,1))),0,0];
height_offset = min(V(:,3));
V(:,3) = V(:,3)-height_offset;
V0(:,3) = V0(:,3)-height_offset;

data.Z = Z;
data.Z0 = Z0;
data.V = V;
data.V0 = V0;
data.F = F;
data.F0 = F0;
data.dataZ = get_data_from_height_field(X,Y,Z,II,Z0);
data.dataZ0 = get_data_from_height_field(X,Y,Z0,II,Z0);
data.dataZ.gaussian = data.dataZ.gaussian(:);
data.dataZ0.gaussian = data.dataZ0.gaussian(:);
data.dataZ.gaussian = data.dataZ.gaussian(J);
data.dataZ0.gaussian = data.dataZ0.gaussian(J0);
data.dataZ.small_sv = data.dataZ.small_sv(:);
data.dataZ0.small_sv = data.dataZ0.small_sv(:);
data.dataZ.small_sv = data.dataZ.small_sv(J);
data.dataZ0.small_sv = data.dataZ0.small_sv(J0);
data.dataZ.ruling_line = data.dataZ.ruling_line(J,:);
data.dataZ0.ruling_line = data.dataZ0.ruling_line(J0,:);
data.dataZ.big_sv = data.dataZ.big_sv(:);
data.dataZ0.big_sv = data.dataZ0.big_sv(:);
data.dataZ.big_sv = data.dataZ.big_sv(J);
data.dataZ0.big_sv = data.dataZ0.big_sv(J0);


save([directory,'/all_data.mat'],'data')
writeOBJ([directory,'/objs/','output.obj'],V,F)
writeOBJ([directory,'/objs/','input.obj'],V0,F0)
if boundary
[data.O,data.L] = ordered_outline(F);
[data.O0,data.L0] = ordered_outline(F0);
end

surface_color = [.8,.8,.8];
view_coords = [0 90];
%%% JUST SURFACE
tsurf(F0,V0,'EdgeColor','none',fphong,fsoft,...
    'FaceVertexCData',repmat(surface_color,size(V0,1),1))
hold on
E0 = edges(F0);
%plot_edges(V0,E0,'-k')
axis equal
axis off
grid off
camproj('orth')
set(gcf,'Color','w');
view(view_coords)
camtarget([0 0 0])
lims = [xlim;ylim;zlim];
camlight
figpng([directory,'/pngs/input_height.png']);
hold off
clf
tsurf(F,V,'EdgeColor','none',fphong,fsoft,...
    'FaceVertexCData',repmat(surface_color,size(V,1),1))
hold on
E = edges(F);
%plot_edges(V,E,'-k')
axis equal
xlim(lims(1,:))
ylim(lims(2,:))
zlim(lims(3,:))
axis off
grid off
camproj('orth')
set(gcf,'Color','w');
view(view_coords)
camtarget([0 0 0])
camlight
figpng([directory,'/pngs/output_height.png']);
hold off
clf

%%% RULING LINES
tsurf(F0,V0,'EdgeColor','none',fphong,fsoft,...
    'FaceVertexCData',repmat(surface_color,size(V0,1),1))
hold on
quiver3(V0(:,1),V0(:,2),V0(:,3),data.dataZ0.ruling_line(:,1),data.dataZ0.ruling_line(:,2),...
    data.dataZ0.ruling_line(:,3),'ShowArrowhead','off')
hold off
axis equal
axis off
grid off
camproj('orth')
set(gcf,'Color','w');
view(view_coords)
camtarget([0 0 0])
lims = [xlim;ylim;zlim];
camlight
figpng([directory,'/pngs/input_rulings.png']);
clf
tsurf(F,V,'EdgeColor','none',fphong,fsoft,...
    'FaceVertexCData',repmat(surface_color,size(V,1),1))
hold on
quiver3(V(:,1),V(:,2),V(:,3),data.dataZ.ruling_line(:,1),data.dataZ.ruling_line(:,2),...
    data.dataZ.ruling_line(:,3),'ShowArrowhead','off')
hold off
axis equal
xlim(lims(1,:))
ylim(lims(2,:))
zlim(lims(3,:))
axis off
grid off
camproj('orth')
set(gcf,'Color','w');
view(view_coords)
camtarget([0 0 0])
camlight
figpng([directory,'/pngs/output_rulings.png']);
clf

%%% 2D IMAGE
twod_image_V0 = tsurf(F0,[V0(:,1),V0(:,2),zeros(size(V0,1),1)],'CData',...
    V0(:,3),'EdgeColor','none');
axis equal
xlim(lims(1,:))
ylim(lims(2,:))
zlim(lims(3,:))
axis off
grid off
camproj('orth')
set(gcf,'Color','w');
view(view_coords)
camtarget([0 0 0])
colormap(gray(255))
clim_bw = caxis;

twodimageV0 = twod_image_V0.CData;
figpng([directory,'/pngs/input_image.png']);
clf
twod_image_V = tsurf(F,[V(:,1),V(:,2),zeros(size(V,1),1)],'CData',...
    V(:,3),'EdgeColor','none');
axis equal
xlim(lims(1,:))
ylim(lims(2,:))
zlim(lims(3,:))
axis off
grid off
camproj('orth')
set(gcf,'Color','w');
view(view_coords)
camtarget([0 0 0])
colormap(gray(255))
caxis(clim_bw)
twodimageV = twod_image_V.CData;
figpng([directory,'/pngs/output_image.png']);

%%% KAPPA2 SCALAR FIELD
tsurf(F0,[V0(:,1),V0(:,2),zeros(size(V0,1),1)],'CData',...
    data.dataZ0.small_sv,'EdgeColor','none')
axis equal
xlim(lims(1,:))
ylim(lims(2,:))
zlim(lims(3,:))
axis off
grid off
camproj('orth')
set(gcf,'Color','w');
caxis([0 100])
view(view_coords)
camtarget([0 0 0])
colormap(cbrewer('Greens',10))
colormap(jet)
figpng([directory,'/pngs/input_small_sv.png']);
clf
tsurf(F,[V(:,1),V(:,2),zeros(size(V,1),1)],'CData',...
    data.dataZ.small_sv,'EdgeColor','none')
axis equal
xlim(lims(1,:))
ylim(lims(2,:))
zlim(lims(3,:))
axis off
grid off
camproj('orth')
set(gcf,'Color','w');
caxis([0 100])
view(view_coords)
camtarget([0 0 0])
colormap(cbrewer('Greens',10))
colormap(jet)
figpng([directory,'/pngs/output_small_sv.png']);
clf

%%% KAPPA1 SCALAR FIELD
tsurf(F0,[V0(:,1),V0(:,2),zeros(size(V0,1),1)],'CData',...
    data.dataZ0.big_sv,'EdgeColor','none')
axis equal
xlim(lims(1,:))
ylim(lims(2,:))
zlim(lims(3,:))
axis off
grid off
camproj('orth')
set(gcf,'Color','w');
caxis([0 100])
view(view_coords)
colormap(cbrewer('Greens',10))
colormap(jet)
figpng([directory,'/pngs/input_big_sv.png']);
clf
tsurf(F,[V(:,1),V(:,2),zeros(size(V,1),1)],'CData',...
    data.dataZ.big_sv,'EdgeColor','none')
axis equal
xlim(lims(1,:))
ylim(lims(2,:))
zlim(lims(3,:))
axis off
grid off
camproj('orth')
set(gcf,'Color','w');
caxis([0 100])
view(view_coords)
camtarget([0 0 0])
colormap(cbrewer('Greens',10))
colormap(jet)
figpng([directory,'/pngs/output_big_sv.png']);
clf

if boundary

hold on
    for l = 1:(numel(data.L)-1)
      plot(V([data.O(data.L(l):(data.L(l+1)-1)) data.O(data.L(l))],1),...
          V([data.O(data.L(l):(data.L(l+1)-1)) data.O(data.L(l))],2),...
          'LineWidth',10,'Color',[228,26,28]./256);
    end
    axis equal
    print('-depsc',[directory,'/pngs/output_loops.eps']);
    figpng([directory,'/pngs/output_loops.png']);
    hold off
    clf
    hold on
    for l = 1:(numel(data.L0)-1)
      plot(V0([data.O0(data.L0(l):(data.L0(l+1)-1)) data.O0(data.L0(l))],1),...
          V0([data.O0(data.L0(l):(data.L0(l+1)-1)) data.O0(data.L0(l))],2),...
          'LineWidth',10,'Color',[228,26,28]./256);
    end
    axis equal
    print('-depsc',[directory,'/pngs/input_loops.eps']);
    figpng([directory,'/pngs/input_loops.png']);
    hold off


end




end

