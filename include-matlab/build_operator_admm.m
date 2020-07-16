function [A,NN] = build_operator_admm(X,Y,Z,II)

m = size(Z,1);
n = size(Z,2);
[NN,BXX,BXY,BYX,BYY,BX,BY,C] = get_connectivity(X,Y,Z,II,'hex');
HH = get_hessians(NN,BXX,BXY,BYX,BYY,BX,BY,C,Z);
[~,SS,~] = get_eigenspaces(HH);
JJ = repmat(NN(:),4,1);
I = [];
for i=1:size(NN,2)
    I = [I;repmat(4*i,length(NN(:,i)),1)];
end
III = [I-3;I-2;I-1;I];
values = [BXX;BXY;BYX;BYY];
A = sparse(III,JJ,values,4*size(NN,2),n*m);
ind = [(1:4:(4*size(NN,2)-3))';...
    (2:4:(4*size(NN,2)-2))';...
    (3:4:(4*size(NN,2)-1))';...
    (4:4:(4*size(NN,2)))'];


end