function R=rotmatrix(A,B)
%RA=B
vv = cross(A,B);
ssc = [0 -vv(3) vv(2); vv(3) 0 -vv(1); -vv(2) vv(1) 0];
R = eye(3) + ssc + ssc^2*(1-dot(A,B))/(norm(vv))^2;
end