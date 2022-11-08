function [writedir] = ProjectionProcessing_pass2(paramfile,manstitchpos)
%% Stitching and generating mean
% % Input: Dataset to process and stitch positions
%% 0.1 Toolboxes
addpath(genpath('./utils'))
%% 0.2 Read param file
% read file, make structure array
fid = fopen(paramfile);
infoStruct = textscan(fid, '%s %s','Delimiter','\t','CommentStyle','//');
infoStruct = cell2struct(infoStruct{2},infoStruct{1},size(infoStruct,1));

% assign variables
samplename = infoStruct.samplename;
tmp = split(infoStruct.heightstep,'-');
hs_range = str2double(tmp{1}):str2double(tmp{2});
nhs = length(hs_range);
infofile = infoStruct.infofile;
rawbasedir = infoStruct.rawdatapath;
projpath = infoStruct.projpath;
projpath2 = infoStruct.projpath2;
stripheight = str2double(infoStruct.stripheight);
h5ImagePath = infoStruct.h5ImagePath;
h5AnglePath = infoStruct.h5AnglePath;

if isfield(infoStruct,'ycrop')
ycrop = cell2mat(textscan(infoStruct.ycrop,'%f%f','Delimiter',','));
ycrop = ycrop(1):ycrop(2);
osy = length(ycrop); % output size y
end
%% 0.3 Load measurement info
% read rings and associated scan names
[nrings,ringnames] = ReadRingsAndNames(infofile,hs_range);

% read sizes, angles, and .par file info
datsize = ReadScanSize([rawbasedir ringnames{1,1} filesep ringnames{1,1} '.nxs'],h5ImagePath);
if isfield(infoStruct,'ycrop')
    ycrop = cell2mat(textscan(infoStruct.ycrop,'%f%f','Delimiter',','));
    ycrop = ycrop(1):ycrop(2);
else
    ycrop = 1:datsize;
end
osy = length(ycrop); % output size y
if osy < datsize(2)
    datsize(2) = osy;
end

[angles,ip180] = ReadAngles([rawbasedir ringnames{1,1} filesep ringnames{1,1} '.nxs'],h5AnglePath);

pixsize_mm = ReadPixelSize([rawbasedir ringnames{1,1} filesep ringnames{1,1} '.par']);
%% 0.4 Sort out filtering options
filtertag = infoStruct.projfiltertype;
filterwidth = str2double(infoStruct.projfilterwidth);
if strcmp(filtertag,'paganin')
    photonEnergy = str2double(infoStruct.photonenergy); % entered into param file as keV
    det_dist_mm = str2double(infoStruct.detdist);
    lambda = lambda_from_E(photonEnergy*1e3); % [m]
    lambda_mm = lambda*1e3;
end
if strcmp(filtertag,'paganin')
    filtfunc = @(proj) paganin_filter_stack(proj,filterwidth,pixsize_mm,lambda_mm,det_dist_mm);
elseif strcmp(filtertag,'gauss')
    filtfunc = @(proj) imgaussfilt(proj,filterwidth);
else
    filtfunc = @(proj) proj;
end
%% 0.5 Set up directories
basedir0 = [projpath samplename filesep];
readdir = [basedir0 'stitched_proj_filtered' filesep];
% directory for writing projections:
basedir = [projpath2 samplename filesep];
if not(isfolder(basedir)); mkdir(basedir); end
writedir = [basedir 'proj' filesep];
if not(isfolder(writedir)); mkdir(writedir); end
%% 0.6 Overlap positions
stitchposy = [0,cumsum(manstitchpos)];
fullheight = ceil(datsize(2)+sum(manstitchpos)); 
%% 1.0 Ring correction, filtering, stitching, saving
t = Tiff([readdir 'proj_uf_h' num2str(1) '_p' num2str(1,'%04d') '.tif'], 'r');
tmp = read(t); close(t);
[sy,sx] = size(tmp);
clear tmp

% tiff settings ()
tagstruct.ImageLength = fullheight; % y
tagstruct.ImageWidth = sx; % x
tagstruct.BitsPerSample = 32; % single precission
tagstruct.RowsPerStrip = stripheight; % strip size (faster loading of strips)
tagstruct.SamplesPerPixel = 1;
tagstruct.Compression = Tiff.Compression.None;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Photometric = Tiff.Photometric.LinearRaw;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

% read ring correction images
rproj = zeros(sy,sx,nhs,'single');
for h = 1:nhs
    th = hs_range(h);
    t = Tiff([readdir 'mproj_h' num2str(th) '.tif'], 'r');
    im = read(t); close(t);
    rproj(:,:,h) = single(im-imgaussfilt(im,50,'Padding','symmetric'));
end
rproj(abs(rproj)>0.1) = 0;

% create blending mask
blendmasky = zeros(fullheight,sx,nhs);
for h = 1:nhs
    if h == 1
        blendmasky(floor(stitchposy(h))+1:floor(stitchposy(h))+datsize(2)-1,:,h) = ones(datsize(2)-1,sx);
        tmp = floor(stitchposy(h+1))+1:floor(stitchposy(h))+datsize(2)-1;
        blendmasky(tmp,:,h) = ones(length(tmp),sx).*linspace(1,0,length(tmp))';
    elseif h == nhs
        blendmasky(floor(stitchposy(h))+2:floor(stitchposy(h))+datsize(2),:,h) = ones(datsize(2)-1,sx);
        tmp = floor(stitchposy(h))+1:floor(stitchposy(h-1))+datsize(2)-1;
        blendmasky(tmp,:,h) = ones(length(tmp),sx).*linspace(0,1,length(tmp))';
    else
        blendmasky(floor(stitchposy(h))+2:floor(stitchposy(h))+datsize(2)-1,:,h) = ones(datsize(2)-2,sx);
        tmp = floor(stitchposy(h+1))+1:floor(stitchposy(h))+datsize(2)-1;
        blendmasky(tmp,:,h) = ones(length(tmp),sx).*linspace(1,0,length(tmp))';
        tmp = floor(stitchposy(h))+1:floor(stitchposy(h-1))+datsize(2)-1;
        blendmasky(tmp,:,h) = ones(length(tmp),sx).*linspace(0,1,length(tmp))';
    end
end
blendmasky = blendmasky./sum(blendmasky,3);
blendmasky = single(blendmasky);

% start filtering, blending
fprintf('Ring correcting, blending, and saving projections...\n'); tic;
parfor p = 1:ip180
    proj = zeros(sy,sx,nhs,'single');
    for h = 1:nhs
        t = Tiff([readdir 'proj_uf_h' num2str(hs_range(h)) '_p' num2str(p,'%04d') '.tif'], 'r');
        proj(:,:,h) = read(t); close(t);
    end
    proj = proj-rproj;
    proj = flipud(proj);
    
    ebump = 10;
    for h = fliplr(1:nhs-1)
        olregi = floor(stitchposy(h+1))+1:floor(stitchposy(h))+datsize(2);
        olregi = datsize(2)-length(olregi)+1:datsize(2);
        olregi = olregi(ebump):olregi(end-ebump);
        olregip1 = (1:length(olregi))+ebump;

        overlapMean = mean(proj(olregi,:,h),'all');
        overlapMeanP1 = mean(proj(olregip1,:,h+1),'all');
        
        %tmp1(h) = overlapMean;
        %tmp2(h) = overlapMeanP1;
        
        proj(:,:,h) = (overlapMeanP1/overlapMean)*proj(:,:,h);
    end
    
    fullproj = zeros(fullheight,sx,nhs);
    for h = 1:nhs
        fullproj(floor(stitchposy(h))+1:floor(stitchposy(h))+datsize(2),:,h) = proj(:,:,h);
    end
    fullproj = sum(fullproj.*blendmasky,3); 
    fullproj = flipud(fullproj);
        
    fullproj = filtfunc(fullproj);
    fullproj = single(fullproj);
    t = Tiff([writedir 'proj_f_' num2str(p,'%04d') '.tif'], 'w');
    t.setTag(tagstruct); t.write(fullproj); t.close();
end
toc
end
%%

