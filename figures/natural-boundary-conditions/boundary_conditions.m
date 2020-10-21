% Replicates the results in Figure 24 completely: will output the exact
% image in the figure

clear;
addpath(genpath('../../'));
rng(4)
hx = 0.05;
hy = hx*sin(pi/3);
d = hx*cos(pi/3);
[X,Y] = meshgrid(-5:hx:5,5:-hy:-5);
m = size(X,1);
n = size(X,2);
X(2:2:m-1,:) = X(2:2:m-1,:)-d;
Z = (0.15.*X.^2)+0.5;
%Z = -sqrt(X.^2+Y.^2) + 5;
Z = Z+.2.*rand(size(X,1),size(X,2));
Z((X(:).^2+Y(:).^2)>(4.7)^2) = 0;
Z0{1} = Z; % Z with noise and full

num_circles = 15;
C = [2,2.5,0.8;...
    -3,-2,0.4;...
    0.2,-0.3,0.6;...
    -3,2,1;...
    2,-3,.8;...
    3,0.1,1;...
    -1,-3,.7];
C(:,3) = C(:,3).*0.2;
for i=1:size(C,1)
    c = (rand(1,2)-.5);
    c = c./normrow(c);
    c = c.*5.*rand(1,1);
    cx = c(1);
    cy = c(2);
    cr = .2.*rand(1,1);
    cx = C(i,1);
    cy = C(i,2);
    cr = C(i,3);
    Z(((X(:)-cx).^2+(Y(:)-cy).^2)<(cr)^2) = 0;
    Z0{2} = Z;
end
surf(X,Y,Z,'EdgeAlpha',0,fsoft,fphong)
Zsol = cell(1,2);
for i=1:numel(Z0)
    [Zsol{i},data] = sparsify_height_field_admm(X,Y,Z0{i},'GetEnergy',false,'UseMex',...
        true,'AggregateNorm',1,'Lambda',100);
end


%save('test_server.mat');
%load('test_server.mat');
for i=1:2
    system(['mkdir bc_cone',num2str(i)]);
    %save_everything(X,Y,Zsol{i},Z0{i},['bc_cone',num2str(i)]);
    [V0{i},F0{i}] = readOBJ(['bc',num2str(i),'/objs/input.obj']);
    [Vsol{i},Fsol{i}] = readOBJ(['bc',num2str(i),'/objs/output.obj']);
end

clf
hold off
clear('t')
off = [0,0,0];
for i=1:numel(Vsol)
    %hold on
    %t{2*(i-1)+1} = tsurf(F0{i},V0{i}+2.*off,falpha(1,0),fsoft,fphong);
    t = tsurf(F0{i},V0{i},falpha(1,0),fsoft,fphong);
    axis equal
    drawnow
    colormap(cbrewer('RdYlBu',10))
    grid off
    axis off
    set(gcf,'Color','w');
    add_isolines(t);
    view([0 90])
    drawnow
    
    figpng(['v0',num2str(i),'.png']);
    %t{2*i} = tsurf(Fsol{i},Vsol{i}+2.*off+[11,0,0],falpha(1,0),fsoft,fphong);
    t = tsurf(Fsol{i},Vsol{i},falpha(1,0),fsoft,fphong);
    axis equal
    drawnow
    colormap(cbrewer('RdYlBu',10))
    grid off
    axis off
    set(gcf,'Color','w');
    add_isolines(t);
    view([0 90])
    drawnow
    figpng(['vsol',num2str(i),'.png']);
    drawnow
    off = off+[11,0,0];
    axis equal
end
colormap(cbrewer('RdYlBu',10))
grid off
axis off
set(gcf,'Color','w');
    

% for i=1:numel(t)
%     add_isolines(t);
% end


