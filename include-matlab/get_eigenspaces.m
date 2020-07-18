function [UU,SS,VV] = get_eigenspaces(HH)
% Singular value decomposition of the Hessian matrices.
    UU = zeros(size(HH));
    SS = UU;
    VV = UU;
    for i=1:size(HH,3)
        [UU(:,:,i),SS(:,:,i),VV(:,:,i)] = svd(HH(:,:,i));
%         if SS(1,1,i)<1e-4
%             SS(1,1,i)=1e-4;
%          end
    end
end
