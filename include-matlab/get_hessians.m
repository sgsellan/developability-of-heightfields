function [HH,normals,ind_term] = get_hessians(NN,BXX,BXY,BYX,BYY,BX,BY,C,Z)
% Given connectivity and information about the grid, build the Hessian matrices.
% This is coded this way so this function is independent of grid type.

HH = zeros(2,2,size(NN,2));
normals = zeros(size(NN,2),3);
ind_term = zeros(size(NN,2),3);

for i=1:size(NN,2)
    z = Z(NN(:,i));
    b = [BXX((size(NN(:,i))*(i-1)+1):(size(NN(:,i))*i))';...
        BXY((size(NN(:,i))*(i-1)+1):(size(NN(:,i))*i))';...
        BYX((size(NN(:,i))*(i-1)+1):(size(NN(:,i))*i))';...
        BYY((size(NN(:,i))*(i-1)+1):(size(NN(:,i))*i))';
        BX((size(NN(:,i))*(i-1)+1):(size(NN(:,i))*i))';...
        BY((size(NN(:,i))*(i-1)+1):(size(NN(:,i))*i))'
        C((size(NN(:,i))*(i-1)+1):(size(NN(:,i))*i))'];
    c = b*z;
    HH(:,:,i) = [c(1),c(2);c(3),c(4)];
    normals(i,:) = [c(5),c(6),-1];
    ind_term(i) = c(7);
end

% Quadric approximation:
% v'*HH(:,:,ind)*v+normals(ind,1:2)*v+ind_term(ind)


end
