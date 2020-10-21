% Replicates the result in Figure 8 by saving different versions
% of a piecewise developable Washington bust for different values
% of our parameter Lambda

addpath(genpath('../../'));
load('../../data/washington.mat'); % load input
omegas = {10,100,1000,10000,100000,1000000,10^7,10^8};

lambda = [];
kappa2 = [];
for i=1:numel(omegas)
    dir_name = ['omega_',num2str(omegas{i})];
%     system(['mkdir ',dir_name]);
     Z = sparsify_height_field_admm(X,Y,Z0,'GetEnergy',false,'UseMex',...
             true,'AggregateNorm',1,...
             'Lambda',omegas{i},'Fill',false,'Jumps',false);
      save_everything(X,Y,Z,Z0,dir_name);
      load([dir_name,'/all_data.mat']);
      kappa2 = [kappa2, median(data.dataZ.small_sv(:))];
      lambda = [lambda,omegas{i}];
      loglog(lambda,kappa2)
      drawnow
end

kappa2 = [kappa2,median(data.dataZ0.small_sv(:))];
lambda = [lambda,10^9];
loglog(lambda,kappa2)
      drawnow
      
      print('-depsc','lambdathing.eps');