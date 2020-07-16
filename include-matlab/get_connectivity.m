function [NN,BXX,BXY,BYX,BYY,BX,BY,C,bad,non_occluded] = get_connectivity(X,Y,Z,II,grid_type,thresh)
m = size(Z,1);
n = size(Z,2);
hx = X(1,2)-X(1,1);
h = hx;
hy = hx;
non_occluded = ones(m*n,1);
bad = [];
if nargin<6
thresh = 40000;
%warning('Changed threshold. 20000')
end

if size(II,2)==1
    II = II'; %row vector
end

switch grid_type
    case 'square'
        NN = [II-1-m;...
            II-1;...
            II-1+m;...
            II-m;...
            II;...
            II+m;...
            II+1-m;...
            II+1;...
            II+1+m];
        A = [hx^2,-hx*hy,hy^2,-hx,hy,1;...
            0,0,hy^2,0,hy,1;...
            hx^2,hx*hy,hy^2,hx,hy,1;...
            hx^2,0,0,-hx,0,1;...
            0,0,0,0,0,1;...
            hx^2,0,0,hx,0,1;...
            hx^2,hx*hy,hy^2,-hx,-hy,1;...
            0,0,hy^2,0,-hy,1;...
            hx^2,-hx*hy,hy^2,hx,-hy,1];
        B = inv(A'*A)*(A');
        BXX = repmat(2.*B(1,:)',length(II),1);
        BXY = repmat(B(2,:)',length(II),1);
        BYX = repmat(B(2,:)',length(II),1);
        BYY = repmat(2.*B(3,:)',length(II),1);
        BX = repmat(B(4,:)',length(II),1);
        BY = repmat(B(5,:)',length(II),1);
        C = repmat(B(6,:)',length(II),1);
        
    case 'hex'
        A = [h^2/4,3*(h^2)/4,-sqrt(3)*(h^2)/4,-0.5*h,sqrt(3)*h/2,1;...
            h^2/4,3*(h^2)/4,sqrt(3)*(h^2)/4,0.5*h,sqrt(3)*h/2,1;...
            h^2,0,0,-h,0,1;...
            0,0,0,0,0,1;...
            h^2,0,0,h,0,1;...
            h^2/4,3*(h^2)/4,sqrt(3)*(h^2)/4,-0.5*h,-sqrt(3)*h/2,1;...
            h^2/4,3*(h^2)/4,-sqrt(3)*(h^2)/4,0.5*h,-sqrt(3)*h/2,1];
        NN = [];
        EE = zeros(6,length(II),6,2);
        
        EE(:,:,1,1) = [II-2;II+m;II-1;II;II-m-1;II+1];
        EE(:,:,2,1) = [II-m;II-2;II;II+m-1;II+m+1;II+2*m-1];
        EE(:,:,3,1) = [II-1-m;II+m-1;II-m;II;II-m+1;II+m+1];
        EE(:,:,4,1) = [II-1;II+2*m-1;II;II+m;II+1;II+2*m+1];
        EE(:,:,5,1) = [II-m+1;II-1;II+1;II;II+2;II+m];
        EE(:,:,6,1) = [II+m-1;II+2*m+1;II;II+m+1;II-m;II+2];
        
        EE(:,:,1,2) = [II-2;II+m;II-m-1;II;II-2*m-1;II-m+1];
        EE(:,:,2,2) = [II-m;II-2;II;II-1;II+1;II+m-1];
        EE(:,:,3,2) = [II-1-2*m;II-1;II-m;II;II-2*m+1;II+1];
        EE(:,:,4,2) = [II-m-1;II+m-1;II;II+m;II-m+1;II+m+1];
        EE(:,:,5,2) = [II-2*m+1;II-m-1;II-m+1;II;II+2;II+m];
        EE(:,:,6,2) = [II-1;II+m+1;II;II+1;II-m;II+2];
        
        a = [hx/2,0];
        b = h*[cos(pi/3),cos(pi/3)];
        pos = [-3*a+b;a+b;-a;a;a-b;3*a-b];
        XX = pos(:,1);
        YY = pos(:,2);
        edge_A = [XX.^2,YY.^2,XX.*YY,XX,YY,ones(size(pos,1),1)];
        edge_B = inv(edge_A'*edge_A)*(edge_A');
        
        
        
        
        for i=1:length(II)
            if mod(mod(II(i),m),2)==1
                NN = [NN,[II(i)-1;...
                    II(i)-1+m;...
                    II(i)-m;...
                    II(i);...
                    II(i)+m;...
                    II(i)+1;...
                    II(i)+1+m]];
%                                 hold off
%                                 surf(X,Y,Z)
%                                 hold on
%                                 plot3(X(II(i)),Y(II(i)),Z(II(i)),'.r','MarkerSize',60)
%                                 plot3(X(NN(:,i)),Y(NN(:,i)),Z(NN(:,i)),'.b','MarkerSize',30)
%                                 view([0 90])
%                                 axis equal
%                                 drawnow
%                                 pause
                if nargout>8
                for j=1:6
                    vals = Z(EE(:,i,j,1));
                    hess = edge_B*vals;
                    S = svd([2*hess(1),hess(3);hess(3),2*hess(2)]);
                    if S(1)>thresh
                        bad = [bad,NN(j,i)];
                        non_occluded(NN(j,i))=0;
                    end
                end
                end
                
            else
                NN = [NN,[II(i)-1-m;...
                    II(i)-1;...
                    II(i)-m;...
                    II(i);...
                    II(i)+m;...
                    II(i)+1-m;...
                    II(i)+1]];
%                                 hold off
%                                 surf(X,Y,Z)
%                                 hold on
%                                 plot3(X(II(i)),Y(II(i)),Z(II(i)),'.r','MarkerSize',60)
%                                 plot3(X(NN(:,i)),Y(NN(:,i)),Z(NN(:,i)),'.b','MarkerSize',30)
%                                 view([0 90])
%                                 axis equal
%                                 drawnow
%                                 pause
                if nargout>8
                for j=1:6
                    vals = Z(EE(:,i,j,2));
                    hess = edge_B*vals;
                    S = svd([2*hess(1),hess(3);hess(3),2*hess(2)]);
                    if S(1)>thresh
                        bad = [bad,NN(j,i)];
                        non_occluded(NN(j,i))=0;
                    end
                end
                end
            end
        end
        B = inv(A'*A)*(A');
        BXX = repmat(2.*B(1,:)',length(II),1);
        BXY = repmat(B(3,:)',length(II),1);
        BYX = repmat(B(3,:)',length(II),1);
        BYY = repmat(2.*B(2,:)',length(II),1);
        BX = repmat(B(4,:)',length(II),1);
        BY = repmat(B(5,:)',length(II),1);
        C = repmat(B(6,:)',length(II),1);
    case 'pinched_hex'
        
        A = [h^2/4,3*(h^2)/4,-sqrt(3)*(h^2)/4,-0.5*h,sqrt(3)*h/2,1;...
            h^2/4,3*(h^2)/4,sqrt(3)*(h^2)/4,0.5*h,sqrt(3)*h/2,1;...
            h^2,0,0,-h,0,1;...
            h^2,0,0,h,0,1;...
            h^2/4,3*(h^2)/4,sqrt(3)*(h^2)/4,-0.5*h,-sqrt(3)*h/2,1;...
            h^2/4,3*(h^2)/4,-sqrt(3)*(h^2)/4,0.5*h,-sqrt(3)*h/2,1];
        NN = [];
        for i=1:length(II)
            if mod(mod(II(i),m),2)==1
                NN = [NN,[II(i)-1;...
                    II(i)-1+m;...
                    II(i)-m;...
                    II(i)+m;...
                    II(i)+1;...
                    II(i)+1+m]];
%                                 hold off
%                                 surf(X,Y,Z)
%                                 hold on
%                                 plot3(X(II(i)),Y(II(i)),Z(II(i)),'.r','MarkerSize',60)
%                                 plot3(X(NN(:,i)),Y(NN(:,i)),Z(NN(:,i)),'.b','MarkerSize',30)
%                                 view([0 90])
%                                 drawnow
%                                 pause
            else
                NN = [NN,[II(i)-1-m;...
                    II(i)-1;...
                    II(i)-m;...
                    II(i)+m;...
                    II(i)+1-m;...
                    II(i)+1]];
                %                 hold off
                %                 surf(X,Y,Z)
                %                 hold on
                %                 plot3(X(II(i)),Y(II(i)),Z(II(i)),'.r','MarkerSize',60)
                %                 plot3(X(NN(:,i)),Y(NN(:,i)),Z(NN(:,i)),'.b','MarkerSize',30)
                %                 view([0 90])
                %                 drawnow
                %                 pause
            end
        end
        B = inv(A'*A)*(A');
        BXX = repmat(2.*B(1,:)',length(II),1);
        BXY = repmat(B(3,:)',length(II),1);
        BYX = repmat(B(3,:)',length(II),1);
        BYY = repmat(2.*B(2,:)',length(II),1);
        
end

% 
% hold off
% surf(X,Y,Z)
% hold on
% plot3(X(bad),Y(bad),Z(bad),'.b','MarkerSize',5)
end