% Running this script will generate all the data in Figure 23
% From left to right, you will find the .obj files rendered in
% that figure at  
% ./bayarea-files/objs/input.obj, 
% ./bayarea-files/objs/output.obj, 
% ./mountain-range-files/objs/input.obj, 
% ./mountain-range-files/objs/output.obj, 
%
% For the specific rendering style, import these objs to 
% ../../render/render-sample.blend and you will exactly replicate 
% the figure

addpath(genpath('../../'));
load('../../data/bayarea.mat') % loading input
Z = sparsify_height_field_admm(X,Y,Z0,'GetEnergy',false,'UseMex',...
    true,'AggregateNorm',1,...
    'Lambda',1000000,'Fill',false); % running method
mkdir('bayarea-files') % saving output (this can take a long time)
save_everything(X,Y,Z,Z0,'bayarea-files')
clear;

load('../../data/range.mat') % loading input
Z = sparsify_height_field_admm(X,Y,Z0,'GetEnergy',false,'UseMex',...
    true,'AggregateNorm',1,...
    'Lambda',1000000,'Fill',false); % running method
mkdir('mountain-range-files') % saving output (this can take a long time)
save_everything(X,Y,Z,Z0,'mountain-range-files')