% This code will generate randomly aligned inputs for our method's 
% results in Fig. 22. Each input and output files will be in subdirectories
% ./rand_*.obj and they can be rendered in Blender as shown in
% ../../render/render-sample.blend


addpath(genpath('../../'))
clear
clf
hx = 0.2;
hy = hx*sin(pi/3);
d = hx*cos(pi/3);
[X,Y] = meshgrid(-5:hx:5,5:-hy:-5);
m = size(X,1);
n = size(X,2);
X(2:2:m-1,:) = X(2:2:m-1,:)-d;
hold off




rng(0) % to ensure same output as in figure (remove this line to try other sets of angles)
offZ = 0;

for ss=1:10
    alpha = rand(1)*2*pi;
    disp(ss)
    disp(alpha*360/(2*pi))
    newX = (1/sqrt(2)).*(X+Y);
    newY = (1/sqrt(2)).*(X-Y);
    V = [newX(:),newY(:)];
    R = [cos(alpha),sin(alpha);-sin(alpha),cos(alpha)];
    V = V*R;
    newX = reshape(V(:,1),m,n);
    newY = reshape(V(:,2),m,n);

Z0 = 12+X.*sin(newY) - Y.*cos(newX);
Z(:,:,1) = -sqrt(10- 0.1.*(newX+newY).^2);
Z(:,:,2) = -10.*ones(size(X,1),size(X,2));
Z(:,:,3) = -sqrt((newX-1).^2+(newY+1).^2);
Z(:,:,4) = sqrt(50 -2.*newY.^2)-9;
Z(:,:,5) = (newX+newY).^2 - 100;
for i=1:size(X,1)
    for j=1:size(X,2)
        Z0(i,j) = max(Z(i,j,:));
    end
end
Z = Z0+6;
Z = Z+.1.*rand(size(X,1),size(X,2));
Z = Z+.05.*sin(5.*newX);
Z((X(:).^2+Y(:).^2)>(4.7)^2) = 0;
     %%%%%%% RUN METHOD %%%%%%%%%%%%
[ZZ1,dataZZ1] = sparsify_height_field_admm(X,Y,Z,'GetEnergy',false,'UseMex',...
    false,'AggregateNorm',1,...
    'Lambda',100,'Fill',false,'Plot',false);

system(['mkdir ','rand',num2str(ss)]);
%save_everything(X,Y,ZZ1,Z,['rand',num2str(ss)]);

II = interior_indeces(Z,1);
[V1,F1] = triangulate_height_field(X,Y,Z,II);
[V2,F2] = triangulate_height_field(X,Y,ZZ1,II);
data = get_data_from_height_field(X,Y,Z,II,Z);
dataZZ = get_data_from_height_field(X,Y,ZZ1,II,Z);
off = [1.4*(max(V1(F1,1))-min(V1(F1,1))),0,0];
%clf

R = [cos(alpha),-sin(alpha),0;sin(alpha),cos(alpha),0;0,0,1];
V1_rot = V1*R;
V2_rot = V2*R;


writeOBJ(['rand',num2str(ss),'_input','.obj'],V1+[0,0,offZ],F1)
writeOBJ(['rand',num2str(ss),'_output','.obj'],V2+off+[0,0,offZ],F2)
%writeOBJ(['rand',num2str(ss),'_output_rotated','.obj'],V2_rot+2.*off+[0,0,offZ],F2)


title(dataZZ1.title_str)
tsurf(F1,V1+[0,0,offZ],'CData',data.gaussian,fphong,fsoft);
hold on
tsurf(F2,V2+off+[0,0,offZ],'CData',dataZZ.gaussian,fphong,fsoft);
tsurf(F2,V2_rot+2.*off+[0,0,offZ],'CData',dataZZ.gaussian,fphong,fsoft);
%[V3,F3]=readOBJ('out.obj');
%t{3} = tsurf(F3,V3+off,'CData',zeros(size(V3,1),1),'EdgeColor',0.5*[1 1 1],fphong);
axis equal
grid off
axis off
colormap(cbrewer('Greens',20))
%caxis([min(min(data.small_sv)),max(max(data.small_sv))]) 
caxis([0,10^8]) 
%colorbar
view([0 67]) 
camproj('persp');
set(gcf,'Color','w');
drawnow
title(dataZZ1.title_str,'FontSize',40)
% t{1}.CData = data.global_energy;
% t{2}.CData = dataZZ.global_energy;
%
%
offZ = offZ+10;
%pause
end

% 
% [V,I] = remove_unreferenced(V2,F2);
% F = I(F);



