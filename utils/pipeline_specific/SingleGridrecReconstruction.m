function [reco] = SingleGridrecReconstruction(name,sino,angles,pixsize_mm,python_path)
%single_gridrec_reconstruction Reconstructs a single sinogram using gridrec
%   This function calls a python script to do the reconstruction
%   Inputs:
%       name - location and name (needs file path!) of files related to this sino/reco
%               sino will be name.h5, reco will be name_reco.tif
%       sino - sinogram dimensions should be (x,theta). COR is middle.
%       angles - list of angles, should be in degrees
%       pixsize_mm - pixel size in mm
%       python_path - path to python interpreter

fprintf('Running gridrec reconstruction...\n')

sino = single(sino);
angles = single(angles);
pixsize_mm = single(pixsize_mm);

if exist([name '.h5'],'file') == 2
    fprintf(['deleting existing temp file ' name '.h5 \n'])
    delete([name '.h5']); 
end 
if exist([name '_reco.h5'],'file') == 2 
    fprintf(['deleting existing temp file ' name '_reco.h5  \n'])
    delete([name '_reco.h5']); 
end 

h5create([name '.h5'],'/sino',size(sino),'Datatype','single')
h5create([name '.h5'],'/angles',size(angles),'Datatype','single')
h5create([name '.h5'],'/pixsize_mm',[1,1],'Datatype','single')
h5write([name '.h5'],'/sino',sino)
h5write([name '.h5'],'/angles',angles)
h5write([name '.h5'],'/pixsize_mm',pixsize_mm)

python_reco_script = 'utils/SingleGridrecReconstruction.py';
CMD=sprintf('%s %s %s',...
    python_path,...
    python_reco_script,...
    [name '.h5']);

[~,result]=system(CMD);

fprintf(result)

t1 = tic;
maxtime = 10;
dt = 0.5;
while ~isfile([name '_reco.h5']) && toc(t1) < maxtime
    pause(dt);
    %elapsedtime = elapsedtime+dt;
end
reco = h5read([name '_reco.h5'],'/reco');
end

