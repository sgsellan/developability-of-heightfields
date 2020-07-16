function [Z,data] = sparsify_height_field_admm(X,Y,Z0,varargin)
% **Find the closest-to-developable surface within an eps of X,Y,Z0**
% Explain how we do it
%
% Input:
%       X,Y are m by n matrices of grid coordinates like those created with
%           MATLAB's meshgrid function
%       Z0 is a m by n matrix containing the Z0 of the input height field
%       Optional:
%
%%%%% PARAMETER ASSIGNING
%%% Default Parameters
II = [];
omega = 10;
maxinneriter = 10000;
get_energy = false;
usemex = true;
aggregatenorm = 1;
opnorm = false;
fill = true;
method = 'admm';
Z00 = [];
jumps = false;
plot = true;
thresh = 400000;
interpolate_b = [];
interpolate_bc = [];
normal_weight = false;
update = false;
%%% End of Default Parameters

params_to_variables = containers.Map({'Interior','Omega',...
    'MaxInnerIter','GetEnergy','UseMex','AggregateNorm',...
    'OperatorNorm','Fill','InitialForGradient','Method','Jumps','Plot','Thresh',...
    'InterpolateB','InterpolateBC','NormalWeight','Update'},...
    {'II','omega','maxinneriter','get_energy',...
    'usemex','aggregatenorm','opnorm','fill','Z00','method','jumps','plot','thresh',...
    'interpolate_b','interpolate_bc','normal_weight','update'});
v = 1;
while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
        assert(v+1<=numel(varargin));
        v = v+1;
        feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
        error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
end


tic

m = size(X,1);
n = size(X,2);
always_II = [];
Z0(Z0==-Inf) = 0;
for j=3:(n-2) % for sub2ind but let's do one thing at a time.
    always_II = [always_II, ((j-1)*m)+(3:(m-2))];
end

remove_ghost = false;
if isempty(II)
    background = unique([find(Z0(:)==0),find(isnan(Z0(:))),find(Z0(:)==-Inf)]);
    II = setdiff(1:m*n,background);
    %II = setdiff(II,interpolate_b);
    remove_ghost = true;
end
II = intersect(II,always_II);
[NN,BXX,BXY,BYX,BYY,BX,BY,C,bad,non_occluded] = ...
    get_connectivity(X,Y,Z0,II,'hex',thresh);
[HH,normals] = get_hessians(NN,BXX,BXY,BYX,BYY,BX,BY,C,Z0);
hess_big_sv = [];
for i=1:size(HH,3)
    [~,s,~] = svd(HH(:,:,i));
    hess_big_sv(i) = s(1,1);
end
background = setdiff(1:(m*n),NN(:))';
ghost = setdiff(NN(:),II);

if jumps

    non_occ = (non_occluded(II));
    remove_ghost = false;
    % non_occ = (hess_big_sv<1000);
    II = II(logical(non_occ));
    NN = NN(:,logical(non_occ));
    background = setdiff(1:(m*n),NN(:))';
    ghost = setdiff(NN(:),II);
    data.II = II;
    if plot
        surf(X,Y,Z0,'EdgeColor','none')
        hold on
        % plot3(X(II),Y(II),Z0(II),'.r','MarkerSize',1)
        plot3(X(background),Y(background),Z0(background),'.g','MarkerSize',1)
        plot3(X(ghost),Y(ghost),Z0(ghost),'.b','MarkerSize',1)
        % drawnow
        hold off
        %  pause
        %%% TESTING SHIT
    end
end
if plot
    surf(X,Y,Z0)
    hold on
    plot3(X(II),Y(II),Z0(II),'.r','MarkerSize',2)
    plot3(X(background),Y(background),Z0(background),'.g','MarkerSize',2)
    plot3(X(ghost),Y(ghost),Z0(ghost),'.b','MarkerSize',2)
    hold off
    drawnow
end


weights = ones(length(II),1);



A = speye(4*length(II),4*length(II));
c = sparse(4*length(II),1);
state.Z = Z0(:);
state.Z0 = Z0(:);
state.U = zeros(size(c,1),size(c,2));
state.rho_prev = nan;
state.rho = .000001;
state.argmin_X_data = [];

Zi = state.Z(1:m*n);
B = build_operator_admm(X,Y,reshape(Zi,m,n),II);

state.X = B*state.Z;
state.argmin_Z_data.rho_prev = nan;
argmin_X = @(Z,U,rho,data) argmin_X_full(Z,U,rho,data,B,II,get_energy,...
    usemex,aggregatenorm,opnorm,weights,NN,BXX,BXY,BYX,BYY,BX,BY,C,update);
argmin_Z = @(X,U,rho,data) argmin_Z_full(X,U,rho,data,B,II,...
    Z0,omega,background,get_energy,weights,interpolate_b,interpolate_bc);


[~,~,state] = admm_changed(argmin_X,argmin_Z,A,-B,c,state,'MaxIter',maxinneriter,...
    'TolAbs',0.0001,'TolRel',0.001);

time = toc;
data.II = II;
data.title_str = ['Lambda: ',num2str(omega),'. Vertices: ',num2str(length(II)),...
    '. Time: ',num2str(time)];
% pause
disp(data.title_str)

if remove_ghost
    state.Z(ghost) = Z0(ghost);
end
%state.Z(ghost) = NaN;

if fill
    Z_gaps = reshape(state.Z,m,n);
    Z_gaps(background) = NaN;
    Z0_gaps = reshape(state.Z0,m,n);
    Z0_gaps(background) = NaN;
    [NN,BXX,BXY,BYX,BYY,BX,BY,C,bad,non_occluded] = ...
        get_connectivity(X,Y,reshape(state.Z,m,n),II,'hex',thresh);
    [HH,normals,ind_term] = get_hessians(NN,BXX,BXY,BYX,BYY,BX,BY,C,reshape(state.Z,m,n));
    interior_verts_X = X(II);
    interior_verts_Y = Y(II);
    interior_verts_Z0 = state.Z0(II);
    interior_verts_initially = [interior_verts_X(:),interior_verts_Y(:),...
        interior_verts_Z0(:)];
    
    for i=1:length(background)
        if state.Z(background(i))>1e-2
            pos = [X(background(i)),Y(background(i)),state.Z(background(i))];
            distances = normrow(interior_verts_initially(:,1:2) - ...
                repmat(pos(:,1:2),length(II),1));
            [~,ind] = min(distances);
            ind = ind(1);
            v = pos-interior_verts_initially(ind,:);
            v = [v(1);v(2)];
            values = .5.*v'*HH(:,:,ind)*v+normals(ind,1:2)*v+ind_term(ind);
            if values>2
                warning('thing')
            end
            state.Z(background(i)) = values;
        else
            state.Z(background(i)) = NaN;
        end
    end
else
    state.Z(background) = state.Z0(background);
end


Z = reshape(state.Z(1:m*n),m,n);
if plot
    surf(X,Y,Z)
    hold on
    plot3(X(II),Y(II),Z(II),'.r','MarkerSize',2)
    plot3(X(background),Y(background),Z0(background),'.g','MarkerSize',2)
    plot3(X(ghost),Y(ghost),Z0(ghost),'.b','MarkerSize',2)
    hold off
    drawnow
end



end