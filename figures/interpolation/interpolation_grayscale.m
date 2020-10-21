%%%%%%%%%%%%%%%%%% DEVELOPABLE INTERPOLATION %%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the closest developable surface given constraints
% Constraints must be given as a *square* .png image (see data/peaks.png)
% Black means "free", every other gray represent constraints.
% There is a 5% padding "margin" which will be ignored.
% Resolution is tied to actual resolution in the image.
%
% For details, see Section 5.3. in "Developability of Heightfields
% via Rank Minimization", ACMSIGGRAPH 2020. Silvia Sellán, 
% Noam Aigerman, Alec Jacobson
%
%%% READING AND PREPARING IMAGE
addpath(genpath('.'))
A = imread('../../data/peaks.png'); % Swap 'peaks.png' for your image
A = A(:,:,1); % we look at R channel only
% Make square
if size(A,1)>size(A,2)
    A = A(1:size(A,2),1:size(A,2));
else
    A = A(1:size(A,1),1:size(A,1));
end
% 
[X,Y,Z0] = image_info_into_hex(double(A)); % convert to hex grid
b = find(Z0(:)>1); % any non-black means constrained
Z0 = Z0./255; % normalize
bc = Z0(b); % gather constraints
Z0 = ones(size(X,1),size(X,2)); % initialization
Z0(b) = bc; % intialize to constraints
Z0(find((X(:)>1.05)+(X(:)<-0.05)+(Y(:)>1.05)+(Y(:)<-0.05))) = 0; % padding
%%%


%%% RUNNING DEVELOPABLE INTERPOLATION
Z = sparsify_height_field_admm(X,Y,Z0,'InterpolateB',b,'InterpolateBC',bc);
%%%




%%% SAVING AND PLOTTING SOLUTION
clf
Z(Z==0) = NaN;
%Z = max(max(Z))-Z;
system('mkdir grayscale');
save_everything(X,Y,Z,Z0,'grayscale',false);
[V,F] = read_triangle_mesh('grayscale/objs/output.obj');
hold off
tsurf(F,V,fsoft,fphong,falpha(1,0));
colormap(cbrewer('RdYlBu',20))
axis equal
axis off
grid off
view([0 26])
camlight
set(gcf,'Color','w');
add_isolines() % if your computer is slow, comment this line
drawnow
%imwrite(myaa({'raw',2}),'grayscale.png') % uncomment to save high-res image




