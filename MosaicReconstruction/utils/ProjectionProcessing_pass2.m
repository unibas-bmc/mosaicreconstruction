function [writedir] = ProjectionProcessing_pass2(paramfile,manstitchpos,...
    translation)
arguments
    paramfile
    manstitchpos
    translation = 0
end
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
if any(translation, 'all')
    manstitchpos = round(manstitchpos - translation(:,1)');
    translation = cumsum(translation(:,2:3));
end
stitchposy = [0,cumsum(manstitchpos)];
fullheight = ceil(datsize(2)+sum(manstitchpos));
%% 1.0 Ring correction, filtering, stitching, saving
tmp = h5read([readdir 'proj_uf_h' num2str(1) '_p' num2str(1,'%04d') ...
    '.h5'], '/proj'); % TODO: optimize e.g. by using h5info
[sy,sx] = size(tmp);
clear tmp

% write angles
h5create([writedir 'angles.h5'], '/angles', size(angles))
h5write([writedir 'angles.h5'], '/angles', angles)

% write angles
if ~isfile([writedir 'angles.h5'])
    h5create([writedir 'angles.h5'], '/angles', size(angles))
    h5write([writedir 'angles.h5'], '/angles', angles)
end

% read ring correction images
rproj = zeros(sy,sx,nhs,'single');
for h = 1:nhs
    th = hs_range(h);
    im = h5read([readdir 'mproj_h' num2str(th) '.h5'], '/proj');
    rproj(:,:,h) = single(im-imgaussfilt(im,50,'Padding','symmetric'));
end
rproj(abs(rproj)>0.1) = 0;

% create blending mask
if nhs == 1
    blendmasky = ones(fullheight,sx,nhs);
else
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
end
blendmasky = single(blendmasky);

% start filtering, blending
fprintf('Ring correcting, blending, and saving projections...\n'); tic;
parfor p = 1:ip180
    proj = zeros(sy,sx,nhs,'single');
    for h = 1:nhs
        proj(:,:,h) = h5read([readdir 'proj_uf_h' num2str(hs_range(h)) ...
            '_p' num2str(p,'%04d') '.h5'], '/proj')
    end
    proj = proj-rproj;
    % translation between height steps in horizontal direction
    if any(translation, 'all')
        alpha = angles(p) * pi / 180.;
        for h = 1:nhs-1
            tx = translation(h,1);
            ty = translation(h,2);
            d_alpha = -ty * cos(alpha) + tx * sin(alpha);
            proj(:,:,h+1) = subpixelshift(proj(:,:,h+1), 0, d_alpha);
        end
    end
    % proj = flipud(proj);
    
    ebump = 1;
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
    
    fullproj = zeros(fullheight,sx,nhs,'single');
    for h = 1:nhs
        fullproj(floor(stitchposy(h))+1:floor(stitchposy(h))+datsize(2),:,h) = proj(:,:,h);
    end
    fullproj = sum(fullproj.*blendmasky,3); 
    
    % fullproj = flipud(fullproj);
        
    fullproj = filtfunc(fullproj);
    fullproj = single(fullproj);
    h5create([writedir 'proj_f_' num2str(p,'%04d') '.h5'], '/proj',...
        size(fullproj), Datatype='single');
    h5write([writedir 'proj_f_' num2str(p,'%04d') '.h5'], '/proj',...
        fullproj)
end
toc
end
%%

