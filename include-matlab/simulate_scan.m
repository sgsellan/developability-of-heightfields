function Z = simulate_scan(Z0,n)

Z = Z0(:);

if nargin==1
    n = 256;
end

minimum = min(min(Z0(Z0~=0)));
maximum = max(max(Z0(Z0~=0)));
bins = linspace(minimum,maximum,n);
for i=1:length(bins)-1
    in_this_bin = find((Z>bins(i)).*(Z<bins(i+1)));
    Z(in_this_bin) = repmat((bins(i)+bins(i+1))/2,length(in_this_bin),1);
end

Z = reshape(Z,size(Z0,1),size(Z0,2));



end

