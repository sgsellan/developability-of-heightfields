% Generates the files shown in Figure 3. They are saved on book/objs/
% with color as .ply files. This can be slow, input scan is very fine


addpath(genpath('../../'));
load('../../data/book.mat'); % load input
n = 1000;
tsurf(F,V,'FaceVertexCData',C(:,1:3),'EdgeColor','none','FaceColor','interp');axis equal;view([-90 -90])
[X,Y,Z,CC] = get_depth_from_viewer(V,F,n,C);
R = CC(:,:,1);
G = CC(:,:,2);
B = CC(:,:,3);
ZZ = sparsify_height_field_admm(X,Y,Z,'GetEnergy',false,'UseMex',...
     true,'AggregateNorm',1,...
     'Lambda',10000000,'Fill',true,'Jumps',false);
 system('mkdir book')
 [J,J0] = save_everything(X,Y,ZZ,Z,'book');
 R = R(J);
 G = G(J);
 B = B(J);
 colors = [R(:),G(:),B(:)];
 
 [V,F] = readOBJ('book/objs/input.obj');
 tsurf(F,V,'FaceVertexCData',colors(:,1:3),'EdgeColor','none','FaceColor','interp');axis equal;view([-90 -90])
 view([-90 -90])
 figpng('book/input.png')
 myWritePLY('book/input_color.ply',V,F,colors(:,1:3));

 
 [V1,F1] = readOBJ('book/objs/output.obj');
 tsurf(F1,V1,'FaceVertexCData',colors(:,1:3),'EdgeColor','none','FaceColor','interp');axis equal;view([-90 -90])
 view([-90 -90])
 figpng('book/output.png')
 myWritePLY('book/output_color.ply',V1,F1,colors(:,1:3));