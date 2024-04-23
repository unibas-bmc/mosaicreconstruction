function [projvol,mprojvol] = LoadProjectionsManualOverlap(paramfile,hs)
%% 0.1 Toolboxes
addpath(genpath('./utils'))
%% 0.2 Read param file
% read file, make structure array
fid = fopen(paramfile);
infoStruct = textscan(fid, '%s %s','Delimiter','\t','CommentStyle','//');
infoStruct = cell2struct(infoStruct{2},infoStruct{1},size(infoStruct,1));

% assign variables
samplename = infoStruct.samplename;
infofile = infoStruct.infofile;
rawbasedir = infoStruct.rawdatapath;
procbasedir = infoStruct.projpath;
verboseMode = logical(str2double(infoStruct.verboseMode)); % shows plots if true
stripheight = str2double(infoStruct.stripheight);
h5ImagePath = infoStruct.h5ImagePath;
h5AnglePath = infoStruct.h5AnglePath;

%% 0.4 Set up directories
procdir = [procbasedir samplename filesep 'tests_hs' num2str(hs,'%02d') filesep];
readdir = [procdir 'proj_manualoverlap' filesep];

%% 0.5 Load measurement info
[nrings,ringnames] = ReadRingsAndNames(infofile,hs);

% read sizes, angles, and .par file info
datsize = ReadScanSize([rawbasedir ringnames{1} filesep ringnames{1} '.nxs'],h5ImagePath);

[angles,ip180] = ReadAngles([rawbasedir ringnames{1} filesep ringnames{1} '.nxs'],h5AnglePath);
ip360 = 2*ip180;

pixsize_mm = ReadPixelSize([rawbasedir ringnames{1} filesep ringnames{1} '.par']);

tmpinfo = imfinfo([readdir 'proj_r' num2str(1) '_' num2str(1,'%04d') '.tif']);
osy = tmpinfo.Height;
%% 1.0 Read projections
projvol = zeros(osy,datsize(2),nrings,ip360,'single');
mprojvol = zeros(osy,datsize(2),nrings,'single');
for i = 1:nrings
    fprintf(['Loading projections for ring ' num2str(i) '/' num2str(nrings) '\n']); 
    tic;
    t = Tiff([readdir 'mrean_proj_r' num2str(i) '.tif'], 'r');
    mprojvol(:,:,i) = read(t);
    parfor p = 1:ip360
        t = Tiff([readdir 'proj_r' num2str(i) '_' num2str(p,'%04d') '.tif'], 'r');
        projvol(:,:,i,p) = read(t);
    end
end
projvol = permute(projvol,[1,2,4,3]);
end

