function data = get_data_from_height_field(X,Y,Z,II,Z0)
% Given a heightfield, obtain principal curvatures, global energy term
% and other quantities for debugging and evaluating outputs in general.

m = size(Z,1);
n = size(Z,2);
hx = abs(X(1,2)-X(1,1));
hy = hx;



Zi = Z(:);
data.gaussian = gaussian_curvature_grid(X,Y,Z);

[NN,BXX,BXY,BYX,BYY,BX,BY,C] = get_connectivity(X,Y,Z,II,'hex');
[Hi,normals] = get_hessians(NN,BXX,BXY,BYX,BYY,BX,BY,C,Z);
        
        
%     data.errors = zeros(m*n,1);
%     data.errors(II) = normrow(abs((A*CCi)-ZZi)');
    
    ruling_line = [];
            Htrue = Hi;
            sv = zeros(m*n,1);
            sv_big = zeros(m*n,1);
       local_energy = zeros(m*n,1);
       global_energy = zeros(m*n,1);
       ruling_line = zeros(size(Htrue,3),2);
       data.ruling_line = zeros(m*n,3);
            for i=1:size(Htrue,3)
                if ~any(isnan(Htrue(:,:,i)))
                [u,sigma,v] = svd(Htrue(:,:,i));
                sv(II(i)) = sigma(2,2);
                sv_big(II(i)) = sigma(1,1);
                [v,eigens] = eigs(Htrue(:,:,i));
                assert(abs(eigens(1,1))>=abs(eigens(2,2)))
                ruling_line(i,:) = v(:,2)';
                local_energy(II(i)) = sigma(2,2)+sigma(1,1);
                end
            end
            global_energy(II) = (Z(II)-Z0(II)).^2;
            data.ruling_line(II,:) = [ruling_line,normals(:,1).*ruling_line(:,1)+...
                normals(:,2).*ruling_line(:,2)];
            data.ruling_line(II,:) = data.ruling_line(II,:)./normrow(data.ruling_line(II,:));
            data.global_energy = reshape(global_energy,m,n);
            data.local_energy = reshape(local_energy,m,n);
            data.big_sv = reshape(sv_big,m,n);
            data.small_sv = reshape(sv,m,n);


end
