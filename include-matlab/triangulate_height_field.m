function [V,F] = triangulate_height_field(X,Y,Z,II)
h = abs(X(1,2)-X(1,1));
if isempty(II)
 II = interior_indeces(Z,h);
end
[NN,BXX,BXY,BYX,BYY] = get_connectivity(X,Y,Z,II,'hex');

Z = Z-min(min(Z(Z~=0)))+.1;

V = [X(:),Y(:),Z(:)];
F = [];

F = [NN(2,:)',NN(1,:)',NN(4,:)';...
        NN(5,:)',NN(2,:)',NN(4,:)';...
        NN(7,:)',NN(5,:)',NN(4,:)';...
        NN(6,:)',NN(7,:)',NN(4,:)';...
        NN(3,:)',NN(6,:)',NN(4,:)';...
        NN(1,:)',NN(3,:)',NN(4,:)'];
    
    %remove zero vertices?
    
    
    
are_non_zero = (V(F(:,1),3)>0).*(V(F(:,2),3)>0).*(V(F(:,3),3)>0);    
F(~are_non_zero,:) = [];


end

