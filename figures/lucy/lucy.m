% Replicates results in Figure 25. Input and output will appear in
% ./lucy-files/objs/input.obj and ./lucy-fines/objs/output.obj
% and they can be render as in ../../render/blender-sample.blend

% We will also compute Gaussian curvatures and output them in a .png image
% called gaussians.png
% The colormap will be red instead of blue, it was
% flipped to blue in Photoshop afterwards for the paper figure

addpath(genpath('../../'));
n = 1000;
omega = 10000000;
load('../../data/lucy.mat');
Z0 = Z;
Z = sparsify_height_field_admm(X,Y,Z0,'GetEnergy',false,'UseMex',...
    true,'AggregateNorm',1,...
    'Lambda',omega,'Fill',true,'Jumps',true);
mkdir('lucy-files') % saving output
save_everything(X,Y,Z,Z0,'lucy-files');

% calculating gaussian curvatures and plotting them
hold off
clf
recalculate = false;
[V,F] = read_triangle_mesh('lucy-files/objs/input.obj');
[V1,F1] = read_triangle_mesh('lucy-files/objs/output.obj');
k = discrete_gaussian_curvature(V,F);
k1 = discrete_gaussian_curvature(V1,F1);
t1 = tsurf(F,V,falpha(1,0),fsoft,fphong,'CData',abs(k));
hold on
t2 = tsurf(F1,V1+[0.8,0,0],fsoft,fphong,falpha(1,0),'CData',abs(k1));
axis equal
view([0 90])

cols = getDistorsionColormap();
colormap(cols)
caxis([0 0.001])
%camlight
grid off;axis off;
set(gcf,'Color','w');
AO1 = apply_ambient_occlusion(t1,'Factor',1);
AO2 = apply_ambient_occlusion(t2,'Factor',1);

clf;
hold off;
t1 = tsurf(F,V,falpha(1,0),fsoft,fphong,'CData',abs(k));
hold on
t2 = tsurf(F1,V1+[0.8,0,0],fsoft,fphong,falpha(1,0),'CData',abs(k1));
axis equal
view([0 90])

cols = getDistorsionColormap();
colormap(cols)
caxis([0 0.001])
grid off;axis off;
set(gcf,'Color','w');

AO1(AO1>0.5)=0.5;
AO2(AO2>0.5)=0.5;
apply_ambient_occlusion(t1,'Factor',1,'AO',AO1);
apply_ambient_occlusion(t2,'Factor',1,'AO',AO2);
figpng('gaussians.png')
