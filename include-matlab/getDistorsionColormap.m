function cols = getDistortionColormap()
	%set to out color map
	CBcols = [0.85 0.85 0.85];%[0.9 0.9 0.9];
	t=(1:64) /64;t=t';
	cols = (1-t)*CBcols + t*[0.7 0 0];
	cols(cols>1) =1;
	caxis([1 10]);
end