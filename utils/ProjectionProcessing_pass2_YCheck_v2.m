function [writedir] = ProjectionProcessing_pass2_YCheck_v2(paramfile,manstitchpos,hs,ycrop1, xcrop)
%% Stitching and generating mean
% % Input: Dataset to process and stitch positions
%% 0.1 Toolboxes
addpath(genpath('./utils'))

%% 0.2 Read param file
hs_range = hs:hs+1;
nhs = length(hs_range);
assert(nhs == 2);

% read file, make structure array
fid = fopen(paramfile);
infoStruct = textscan(fid, '%s %s','Delimiter','\t','CommentStyle','//');
infoStruct = cell2struct(infoStruct{2},infoStruct{1},size(infoStruct,1));

% assign variables
samplename = infoStruct.samplename;
tmp = split(infoStruct.heightstep,'-');
infofile = infoStruct.infofile;
rawbasedir = infoStruct.rawdatapath;
projpath = infoStruct.projpath;
stripheight = str2double(infoStruct.stripheight);
h5ImagePath = infoStruct.h5ImagePath;
h5AnglePath = infoStruct.h5AnglePath;

%% 0.3 Load measurement info and sort out filtering options
% read rings and associated scan names
[nrings,ringnames] = ReadRingsAndNames(infofile,hs_range);

% read sizes, angles, and .par file info
datsize = ReadScanSize([rawbasedir ringnames{1,1} filesep ringnames{1,1} '.nxs'],h5ImagePath);
[angles,ip180] = ReadAngles([rawbasedir ringnames{1,1} filesep ringnames{1,1} '.nxs'],h5AnglePath);

pixsize_mm = ReadPixelSize([rawbasedir ringnames{1,1} filesep ringnames{1,1} '.par']);

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
%% 0.4 Set up directories
basedir0 = [projpath samplename filesep];
readdir = [basedir0 'stitched_proj_filtered' filesep];
% directory for writing projections:
basedir = [projpath samplename filesep];
if not(isfolder(basedir)); mkdir(basedir); end
writedir = [basedir 'proj_y_test' filesep];
if not(isfolder(writedir)); mkdir(writedir); end
%% 0.6 Overlap positions
stitchposy = [0,cumsum(manstitchpos)];
fullheight = ceil(datsize(2)+sum(manstitchpos));

ycrop2 = ycrop1 - round(stitchposy(2));
[strips1, ycrop1_strip] = findTiffStrips(ycrop1, stripheight);
[strips2, ycrop2_strip] = findTiffStrips(ycrop2, stripheight);
%% 1.0 Ring correction, filtering, stitching, saving
t = Tiff([readdir 'proj_uf_h' num2str(hs) '_p' num2str(1,'%04d') '.tif'], 'r');
tmp = read(t); close(t);
[~,sx] = size(tmp);
xcrop = xcrop+1:sx-xcrop;

% tiff settings ()
osy = length(ycrop1);
tagstruct.ImageLength = length(ycrop1); % y
tagstruct.ImageWidth = length(xcrop); % x
tagstruct.BitsPerSample = 32; % single precision
tagstruct.RowsPerStrip = stripheight; % strip size (faster loading of strips)
tagstruct.SamplesPerPixel = 1;
tagstruct.Compression = Tiff.Compression.None;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

% read ring correction images
rproj = zeros(length(ycrop1),length(xcrop),nhs,'single');
for h = 1:nhs
    th = hs_range(h);
    t = Tiff([readdir 'rproj_h' num2str(th) '.tif'], 'r');
    im = read(t); close(t);
    if h == 1
        rproj(:,:,h) = im(ycrop1, xcrop);
    elseif h == 2
        rproj(:,:,h) = im(ycrop2, xcrop);
    end
end
rproj(abs(rproj)>0.1) = 0;

% % create blending mask
% blendmasky = zeros(fullheight,sx,nhs);
% for h = 1:nhs
%     if h == 1
%         blendmasky(floor(stitchposy(h))+1:floor(stitchposy(h))+datsize(2)-1,:,h) = ones(datsize(2)-1,sx);
%         tmp = floor(stitchposy(h+1))+1:floor(stitchposy(h))+datsize(2)-1;
%         blendmasky(tmp,:,h) = ones(length(tmp),sx).*linspace(1,0,length(tmp))';
%     elseif h == nhs
%         blendmasky(floor(stitchposy(h))+2:floor(stitchposy(h))+datsize(2),:,h) = ones(datsize(2)-1,sx);
%         tmp = floor(stitchposy(h))+1:floor(stitchposy(h-1))+datsize(2)-1;
%         blendmasky(tmp,:,h) = ones(length(tmp),sx).*linspace(0,1,length(tmp))';
%     else
%         blendmasky(floor(stitchposy(h))+2:floor(stitchposy(h))+datsize(2)-1,:,h) = ones(datsize(2)-2,sx);
%         tmp = floor(stitchposy(h+1))+1:floor(stitchposy(h))+datsize(2)-1;
%         blendmasky(tmp,:,h) = ones(length(tmp),sx).*linspace(1,0,length(tmp))';
%         tmp = floor(stitchposy(h))+1:floor(stitchposy(h-1))+datsize(2)-1;
%         blendmasky(tmp,:,h) = ones(length(tmp),sx).*linspace(0,1,length(tmp))';
%     end
% end
% blendmasky = blendmasky./sum(blendmasky,3);
% blendmasky = single(blendmasky);

% start filtering, blending
fprintf('Ring correcting, blending, and saving projections...\n'); tic;
parfor p = 1:ip180
    fprintf('Projection %d\n', p);
    proj = zeros(length(ycrop1), length(xcrop), 2, 'single')
    proj_tmp = zeros(length(strips1) * stripheight,sx,'single');
    
    t = Tiff([readdir 'proj_uf_h' num2str(hs_range(1)) '_p' ...
        num2str(p,'%04d') '.tif'], 'r');
    for s = 1:length(strips1)
        stripNumber = strips1(s);
        proj_tmp((s-1)*stripheight+1:s*stripheight,:) = ...
            readEncodedStrip(t, stripNumber);
    end
    proj(:,:,1) = proj_tmp(ycrop1_strip, xcrop)
    close(t);
    
    t = Tiff([readdir 'proj_uf_h' num2str(hs_range(2)) '_p' ...
        num2str(p,'%04d') '.tif'], 'r');
    for s = 1:length(strips2)
        stripNumber = strips2(s);
        proj_tmp((s-1)*stripheight+1:s*stripheight,:) = ...
            readEncodedStrip(t, stripNumber);
    end
    proj(:,:,2) = proj_tmp(ycrop2_strip, xcrop)
    close(t);
    
    proj = proj-rproj;
    
    overlapMean = mean(proj, [1, 2]);
    proj(:,:,2) = overlapMean(1) / overlapMean(2) * proj(:,:,2);
    
    proj = single(proj);
    proj = filtfunc(proj);
    t = Tiff([writedir 'proj1_f_' num2str(p,'%04d') '.tif'], 'w');
    t.setTag(tagstruct); t.write(proj(:,:,1)); t.close();
    t = Tiff([writedir 'proj2_f_' num2str(p,'%04d') '.tif'], 'w');
    t.setTag(tagstruct); t.write(proj(:,:,2)); t.close();
end
toc
end
%%

