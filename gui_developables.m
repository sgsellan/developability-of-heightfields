function gui_developables(filename,grid_size)



if nargin==1
    grid_size = 200;
end

[filepath,name,ext] = fileparts(filename);
datestring = datestr(now,'mmmm-dd-HH-MM-SS');
dir_name = [name,'-',datestring];
objs_dir_name = [dir_name,'/objs'];
pngs_dir_name = [dir_name,'/pngs'];
mkdir_call = ['mkdir ',dir_name];
system(mkdir_call);



switch ext
    case '.obj'
        [V,F]=readOBJ(filename);
    case '.ply'
        [V,F]=readPLY(filename);
    case '.off'
        [V,F]=readOFF(filename);
    case '.stl'
        [V,F]=readSTL(filename);
end
V=V-min(V);
V=V./max(max(V));
writeOBJ([dir_name,'/',name,'_mesh.obj'],V,F);

%%% EXPLANATION
disp('Welcome to the "Developability of Heightfields via Rank Minimization" Graphic Interface')
disp('(1) You can rotate the figure to pick your preferred viewing angle')
disp('(2) Slide the slider to pick a Lambda (more to the right means more faithful to input)')
disp('(3) Choose whether to detect occlusions or not by toggling the button')
disp('(4) Click RUN to run our method')
%%%

%%% Set-up view
tsurf(F,V)
axis equal
camproj('orth')
%%%

%%% Set up buttons from GUI
GetDepthButton = uicontrol('Style', 'pushbutton','Callback',@getdepthbuttonpress,'String','Just get height',...
    'Position',[20    50    120    20]);
RunButton = uicontrol('Style', 'pushbutton','Callback',@runbuttonpress,'String','RUN',...
    'Position',[20    20    120    20]);
LambdaSlider = uicontrol('Style','slider','Value',5,'Min',1,'Max',10,'String','log(\lambda)',...
    'Position',[20    120    120    20]);
jumpstoggle = uicontrol('Style','togglebutton','Value',0,'Min',0,'Max',1,'String','Detect occlusions',...
    'Position',[20    80    120    20]);
%%%

%%% Functions defining what happens when a button is pressed
    function getdepthbuttonpress(src,e)
        [X,Y,Z0] = get_depth_from_viewer(V,F,grid_size);
        surf(X,Y,Z0)
    end
    function runbuttonpress(src,e)
        lambda = 10^(LambdaSlider.Value);
        [X,Y,Z0] = get_depth_from_viewer(V,F,grid_size);
        surf(X,Y,Z0)
        disp('Running ADMM...')
        Z = sparsify_height_field_admm(X,Y,Z0,'GetEnergy',false,'UseMex',...
            true,'AggregateNorm',1,'Omega',lambda,'Fill',true,'Jumps',jumpstoggle.Value);
        disp('Done! Saving data...')
        save_everything(X,Y,Z,Z0,dir_name,false);
    end
%%%
end
