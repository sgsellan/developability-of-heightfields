function [II,Zcorr,ghost] = interior_indeces(Z,h)

ZZ = Z(:);
ZZcorr = Z(:);
m = size(Z,1);
n = size(Z,2);
II =[]; % These three lines build indeces. Should swap all this
for j=2:(n-1) % for sub2ind but let's do one thing at a time.
    II = [II, ((j-1)*m)+(2:(m-1))];
end
interior = (ZZ(II)~=0).*(ZZ(II-1)~=0).*(ZZ(II+1)~=0).*(ZZ(II-m)~=0).*...
    (ZZ(II-m-1)~=0).*(ZZ(II-m+1)~=0).*(ZZ(II+m)~=0).*...
    (ZZ(II+m-1)~=0).*(ZZ(II+m+1)~=0);
%ZZcorr(II(~logical(interior'))) = 0;
II = II(logical(interior'));
all = 1:(m*n);
is_in_stencil = zeros(n*m,1);
is_in_stencil(II) = 1;
is_in_stencil(II+1) = 1;
is_in_stencil(II-1) = 1;
is_in_stencil(II+m) = 1;
is_in_stencil(II-m) = 1;
is_in_stencil(II+1+m) = 1;
is_in_stencil(II+1-m) = 1;
is_in_stencil(II-1+m) = 1;
is_in_stencil(II-1-m) = 1;
background = find(is_in_stencil == 0);
ZZcorr(background) = 0;

not_background = setdiff(1:(m*n),background')';
nb_stencil = zeros(n*m,1);
nb_stencil(not_background) = 1;
nb_stencil(not_background+1) = 1;
nb_stencil(not_background-1) = 1;
nb_stencil(not_background+m) = 1;
nb_stencil(not_background-m) = 1;
nb_stencil(not_background+1+m) = 1;
nb_stencil(not_background+1-m) = 1;
nb_stencil(not_background-1+m) = 1;
nb_stencil(not_background-1-m) = 1;
nb_stencil = find(nb_stencil>0);
ghost = setdiff(nb_stencil,not_background');


% neigh = [ZZ(II-1-m)';...
%         ZZ(II-1)';...
%         ZZ(II-1+m)';...
%         ZZ(II-m)';...
%         ZZ(II)';...
%         ZZ(II+m)';...
%         ZZ(II+1-m)';...
%         ZZ(II+1)';...
%         ZZ(II+1+m)'];
%     
% non_occ = abs(max(neigh)-min(neigh))<7*h;
% II = II(logical(non_occ));


Zcorr = reshape(ZZcorr,m,n);

    
    
end
