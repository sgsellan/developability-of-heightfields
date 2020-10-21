function myWritePLY(filename,V,F,C,colorMapName)
if nargin  < 5
    colorMapName = 'default';
end
mode = 'ascii';
if(isempty(C)) % no color information
    fid = fopen(filename,'wt');
    fprintf(fid,'ply\nformat %s 1.0\n',mode);
    fprintf(fid,'element vertex %d\n',size(V,1));
    fprintf(fid,'property double x\n');
    fprintf(fid,'property double y\n');
    fprintf(fid,'property double z\n');
    fprintf(fid,'element face %d\n',size(F,1));
    fprintf(fid,'property list int int vertex_indices\n');
    fprintf(fid,'end_header\n');
    FF = [size(F,2)*ones(size(F,1),1) F-1];
    fprintf(fid,'%0.17f %0.17f %0.17f %d %d %d\n',VC');
    format = [repmat('%d ',1,size(FF,2)) '\n'];
    fprintf(fid,format,FF');
else % has color
    fid = fopen(filename,'wt');
    fprintf(fid,'ply\nformat %s 1.0\n',mode);
    fprintf(fid,'element vertex %d\n',size(V,1));
    fprintf(fid,'property double x\n');
    fprintf(fid,'property double y\n');
    fprintf(fid,'property double z\n');
    fprintf(fid,'property uchar red\n');
    fprintf(fid,'property uchar green\n');
    fprintf(fid,'property uchar blue\n');
    fprintf(fid,'property uchar alpha\n');
    fprintf(fid,'element face %d\n',size(F,1));
    fprintf(fid,'property list int int vertex_indices\n');
    fprintf(fid,'end_header\n');
    FF = [size(F,2)*ones(size(F,1),1) F-1];
    
    if (size(C,2) == 1) % scalar value function
        CMap = myColorMap(colorMapName);
        cmin = min(C(:));
        cmax = max(C(:));
        m = size(CMap,1);
        idx = floor((C-cmin)/(cmax-cmin)*m)+1; 
        RGB = ind2rgb(idx,CMap);
        RGB = reshape(RGB, size(RGB,1), size(RGB,3));
        RGB = round(RGB * 255);
        RGBA = [RGB, 255*ones(size(RGB,1),1)];   
    elseif (size(C,2) == 3) % RGB color
        if (max(C(:),1) <= 1)
            RGB = round(C * 255);
        else
            RGB = round(C);
        end
        RGBA = [RGB, 255*ones(size(RGB,1),1)];
    elseif (size(C,2) == 4) % RGBA
        if (max(C(:),1) <= 1)
            RGBA = round(C * 255);
        else
            RGBA = round(C);
        end
    end
    VC = [V, RGBA];
    fprintf(fid,'%0.17f %0.17f %0.17f %d %d %d %d\n',VC');
    format = [repmat('%d ',1,size(FF,2)) '\n'];
    fprintf(fid,format,FF');
end
fclose(fid);