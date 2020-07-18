function [Z,data_2] = argmin_Z_full(X,U,rho,data,A,II,Z0,omega,background,...
    get_energy,weights,interpolate_b,interpolate_bc)
% Solves Z update in ADMM, which amounts to
%               min lambda*||Z-Z0||^2 + (rho/2)||U-A*Z-X||^2
%
% The left hand side in the resulting linear system of equations
% is precomputed, and only re-precomputed if rho changed.

quad_error = false;

% II = 1:m*n;
nm2 = size(A,2);
nm = nm2 - length(II);

% m = size(A,2);
% Z = sdpvar(m^2,1,'full');
% objective = omega*norm(Z-Z0,1)+(rho/2)*(norm(X-(A*Z)+U)^2);
% constraints = [Z(background) == 0];
% options = sdpsettings('verbose',0,'debug',1,'solver','mosek',...
%     'cachesolvers',1);
% sol = optimize(constraints,objective,options);
% Z = value(Z);
L1 = false;
if L1
    
data_2 = data;

C = X+U;
% nm2 = size(A,2);
% nm = nm2 - length(II);
% ss = nm2+nm;
% x0 = [Z0(:);zeros(length(II),1);zeros(length(Z0(:)),1)];
% 
% Q = sparse(ss,ss);
% Q(1:nm2,1:nm2) = rho.*A'*A;
% L = sparse(ss,1);
% L(1:nm2) = -rho.*C'*A;
% L(nm2+1:end) = omega;
% backIdiag = zeros(nm2,1);
% backIdiag(background) =1;
% backI = sparse(1:nm2,1:nm2,backIdiag);
% constraints = [-speye(nm2,nm2),-speye(nm2,nm2);...
%     -speye(nm2,nm2),speye(nm2,nm2)];
% rh_constraints = zeros(size(constraints,1),1);
% eq_constraints = [backI,sparse(nm,nm)];
% rh_eq_constraints = zeros(size(eq_constraints,1),1);
% nm2 = size(A,2);
% nm = nm2-length(II);
% ss = nm2+nm;
% Q = sparse(ss,ss);
% Q(1:nm2,1:nm2) = rho.*A'*A; %quadprog adds the .5
% L = sparse(ss,1);
% L(1:nm2) = -rho.*C'*A;
% L(nm2+1:end) = omega;
% eq_constraints = sparse(1:length(background),background,...
%     1,length(background),ss);
% rh_eq_constraints = zeros(size(eq_constraints,1),1);
% just_nm_big = sparse(1:nm,1:nm,1,nm,nm2);
% just_nm_smol = sparse(1:nm,1:nm,1,nm,nm);
% ineq_constraints = [-just_nm_big,-just_nm_smol;...
%     just_nm_big,-just_nm_smol];
% ineq_rhs = [-Z0(:);Z0(:)];
% 
% x0 = [Z0(:);zeros(length(II),1);zeros(length(Z0(:)),1)];
% 
% disp('calling quadprog')
% [uu,fval,exitflag,output,lambda]=quadprog(Q,L,ineq_constraints,ineq_rhs,eq_constraints,rh_eq_constraints,...
%     -inf,inf,x0);
% 
% 
% Z = uu(1:nm2);

eps = 1e-9;


H = rho.*A'*A;
H = H+eps.*speye(size(H));
Hinv = inv(H);
ZZ0 = [Z0(:);zeros(length(II),1)];
u = zeros(size(ZZ0,1),1);
d = C-A*ZZ0;

for s=1:5
    G = rho.*u'*A'*A-rho.*d'*A+omega.*sign(u)';
    u = Hinv*(-G');
   % assert(norm(H*u+G',inf)<1e-1);
end
    

Z = u;






else
%%$$ HERE ONWARDS
nm2 = size(A,2);
nm = nm2 - length(II);
test = false;


ZZ0 = Z0(:);
nm = size(A,2);
C = X+U;
eps = 1e-2;
omega = omega.*weights;
if isempty(interpolate_b)
I = sparse(II,II,omega,nm,nm) + eps.*sparse(1:nm,1:nm,1,nm,nm);
else
    I = eps.*sparse(1:nm,1:nm,1,nm,nm);
end
%I = omega.*sparse(1:nm,1:nm,1,nm,nm);
%I = omega.*sparse(II,II,1,nm,nm);
if rho == data.rho_prev
    data_2.Q = data.Q;
    if quad_error
        data_2.Q = data_2.Q+10^(8).*E'*E;
    end
else
    if test
        data_2.Q = omega.*A'*A+.5.*rho.*A'*A;
        disp('DOING HESSIAN FIDELITY')
    else
    data_2.Q = I+.5.*rho.*A'*A;
    end
end
if test
    L = -2.*((omega.*ZZ0'*A'*A)+(.5.*rho.*C'*A));
else
L = -2.*((ZZ0'*I)+(.5.*rho.*C'*A));
end


if isempty(interpolate_b)
    b = background;
    bc = ZZ0(background);
else
    
    b_repeated = [background;interpolate_b];
    bc_repeated = [ZZ0(background);interpolate_bc];
    [b,IA,IC] = unique(b_repeated);
    bc = bc_repeated(IA);
    
    
end

if rho == data.rho_prev
    [uu,data_2.F] = min_quad_with_fixed(data_2.Q,L',b,bc,...
        [],[],data.F);
    % disp(['thing'])
else
    [uu,data_2.F] = min_quad_with_fixed(data_2.Q,L',b,bc);
end

data_2.rho_prev = rho;
if get_energy
    ind_term = ZZ0'*I*ZZ0+.5.*rho*C'*C;
    data_2.energy = uu'*data_2.Q*uu+L*uu+ind_term;
end

Z = uu;

%disp(max(abs(A*(Z-ZZ0))))
end

end

