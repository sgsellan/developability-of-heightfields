function K = gaussian_curvature_grid(X,Y,Z)
% Gaussian curvature on a square grid calculated via angle defect.

    K = zeros(size(X,1)*size(X,2),1);

    m = size(X,1);
    n = size(X,2);   
    II =[];
    for j=2:(n-1)
        II = [II, ((j-1)*m)+(2:(m-1))];
    end

    VV = [X(:),Y(:),Z(:)];
    
    ax1 = VV(II,:)-VV(II-m,:);
    ax2 = VV(II,:)-VV(II+m,:);
    ay1 = VV(II,:)-VV(II-1,:);
    ay2 = VV(II,:)-VV(II+1,:);
    
    ang1 = acos(dot(ax1,ay1,2)./(normrow(ax1).*normrow(ay1)));
    ang2 = acos(dot(ax1,ay2,2)./(normrow(ax1).*normrow(ay2)));
    ang3 = acos(dot(ax2,ay2,2)./(normrow(ax2).*normrow(ay2)));
    ang4 = acos(dot(ax2,ay1,2)./(normrow(ax2).*normrow(ay1)));
    
    KK = (2*pi) - (ang1+ang2+ang3+ang4);
    
    K(II) = KK;
    
    K = reshape(K,m,n);
end
