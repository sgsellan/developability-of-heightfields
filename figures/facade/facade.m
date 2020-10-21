% Replicates the result in Figure 6 in .obj files which will be
% generated in facade/objs/input.obj and facade/objs/output.obj

addpath(genpath('../../'));
load('../../data/facade_rotated.mat') % Load input
Z0 = Z;
Z = sparsify_height_field_admm(X,Y,Z0,'GetEnergy',false,'UseMex',...
    true,'AggregateNorm',1,...
    'Lambda',1000000,'Fill',false); % run method
system('mkdir facade') % save output
save_everything(X,Y,Z,Z0,'facade')
