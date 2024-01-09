function [writedir] = ProjectionProcessing_pass1(paramfile,manstitchpos)
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
stripheight = str2double(infoStruct.stripheight);
h5ImagePath = infoStruct.h5ImagePath;
h5AnglePath = infoStruct.h5AnglePath;
imagefmt = infoStruct.imagefmt;

if isfield(infoStruct,'ycrop')
ycrop = cell2mat(textscan(infoStruct.ycrop,'%f%f','Delimiter',','));
ycrop = ycrop(1):ycrop(2);
osy = length(ycrop); % output size y
end

if isfield(infoStruct,'zingerCorrection')
    zingerCorrection = str2double(infoStruct.zingerCorrection);
else
    zingerCorrection = 1;
end
if isfield(infoStruct,'gainoffsetCorrection')
    gainoffsetCorrection = str2double(infoStruct.gainoffsetCorrection);
    load(infoStruct.gainoffsetPath);
    
else
    gainoffsetCorrection = 0;
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
basedir = [projpath samplename filesep];
if not(isfolder([basedir])); mkdir([basedir]); end
% directory for writing projections:
writedir = [basedir 'stitched_proj_filtered' filesep];
if not(isfolder(writedir)); mkdir(writedir); end
%% 0.6 Overlap positions
%% If manstitchpos has size [1, nrings], take the same stitching positions
%% for all height steps, otherwise use manstitchpos(h,:) for height step h.
% tmpdir = [procdir 'parameters' filesep];
% load([tmpdir 'stitching_parameters.mat'],'ol_final_subpix','corpix_subpix','fullwidth','stitchpos','halfwidth')
corpix_subpix_l = manstitchpos(:,1);
ol_final_subpix_l = manstitchpos(:,2:end);
halfwidth_l = ceil(datsize(2)+sum(ol_final_subpix_l,2));
fullwidth_l = halfwidth_l*2-ceil(corpix_subpix_l);
fullwidth_min = min(fullwidth_l);
%% 1.0 Loop over heights
for h = 1:nhs
if size(corpix_subpix_l, 1) == 1
    corpix_subpix = corpix_subpix_l;
    ol_final_subpix = ol_final_subpix_l;
    halfwidth = halfwidth_l;
    fullwidth = fullwidth_l;
else
    corpix_subpix = corpix_subpix_l(h,:);
    ol_final_subpix = ol_final_subpix_l(h,:);
    halfwidth = halfwidth_l(h,:);
    fullwidth = fullwidth_l(h,:);
end

stitchpos = [0,cumsum(ol_final_subpix)];
% worst case we crop asymetrically by 1 pixel
crop_left = ceil((fullwidth - fullwidth_min) / 2);
crop_right = fullwidth - fullwidth_min - crop_left;
th = hs_range(h);
fprintf('Working on height %d (%d/%d)\n',th,h,nhs)
%% 1.1 Flats/darks loading
% using mean flat and median dark
mflats = zeros(datsize(1),datsize(2),nrings);
mddarks = zeros(datsize(1),datsize(2),nrings);
for i = 1:nrings
    tmp = rot90(h5read([rawbasedir ringnames{h,i} filesep 'pre_ref.nxs'],'/ref_0'));
    tmp1 = rot90(h5read([rawbasedir ringnames{h,i} filesep 'post_ref.nxs'],'/ref_0'));
    flats = double(cat(3,tmp,tmp1)); clear tmp tmp1
    mflats(:,:,i) = mean(flats,3);
    darks = rot90(h5read([rawbasedir ringnames{h,i} filesep 'post_dark.nxs'],'/dark_0'));
    mddarks(:,:,i) = double(median(darks,3));
end
clear darks flats tmp1 tmp
%% 1.2 Flat/dark correction, gain/offset/zinger, stitching, saving
% tiff settings ()
tagstruct.ImageLength = osy; % y
tagstruct.ImageWidth = fullwidth_min; % x
tagstruct.BitsPerSample = 32; % single precission
tagstruct.RowsPerStrip = stripheight; % strip size (faster loading of strips)
tagstruct.SamplesPerPixel = 1;
tagstruct.Compression = Tiff.Compression.None;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Photometric = Tiff.Photometric.LinearRaw;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

% build a mask for the 0 degrees projection
halfmask = zeros(osy,halfwidth,nrings);
for r = 1:nrings
    if r == 1
        halfmask(:,floor(stitchpos(r))+1:floor(stitchpos(r))+datsize(2)-1,r) = ones(osy,datsize(2)-1);
        tmp = floor(stitchpos(r+1))+1:floor(stitchpos(r))+datsize(2)-1;
        halfmask(:,tmp,r) = ones(osy,length(tmp)).*linspace(1,0,length(tmp));
    elseif r == nrings
        halfmask(:,floor(stitchpos(r))+2:floor(stitchpos(r))+datsize(2),r) = ones(osy,datsize(2)-1);
        tmp = floor(stitchpos(r))+1:floor(stitchpos(r-1))+datsize(2)-1;
        halfmask(:,tmp,r) = ones(osy,length(tmp)).*linspace(0,1,length(tmp));
    else
        halfmask(:,floor(stitchpos(r))+2:floor(stitchpos(r))+datsize(2)-1,r) = ones(osy,datsize(2)-2);
        tmp = floor(stitchpos(r+1))+1:floor(stitchpos(r))+datsize(2)-1;
        halfmask(:,tmp,r) = ones(osy,length(tmp)).*linspace(1,0,length(tmp));
        tmp = floor(stitchpos(r))+1:floor(stitchpos(r-1))+datsize(2)-1;
        halfmask(:,tmp,r) = ones(osy,length(tmp)).*linspace(0,1,length(tmp));
    
    end
end
% build a mask for the 0 plus 180 degrees projection
fullmask = zeros(osy,fullwidth,2);
fullmask(:,1:halfwidth-1,1) = ones(osy,halfwidth-1);
fullmask(:,halfwidth-ceil(corpix_subpix)+2:end,2) = ones(osy,halfwidth-1);
tmp = halfwidth-ceil(corpix_subpix)+2:halfwidth;
fullmask(:,halfwidth-ceil(corpix_subpix)+2:halfwidth,1) = ones(osy,length(tmp)).*linspace(1,0,length(tmp));
fullmask(:,halfwidth-ceil(corpix_subpix)+2:halfwidth,2) = ones(osy,length(tmp)).*linspace(0,1,length(tmp));

% pre-crop flats/darks to avoid repeating crop ip180 times
mddarks_crop = mddarks(ycrop,:,:);
mflats_crop = mflats(ycrop,:,:)-mddarks_crop; % avoid repeating subtraction

% crop gains and offsets
if gainoffsetCorrection == 1
    gains = gains(ycrop,:);
    offsets = offsets(ycrop,:);
else
    gains = ones(osy,datsize(2));
    offsets = zeros(osy,datsize(2));
end

% start filtering, blending, stitching
fprintf('Stitching and saving projections...\n'); tic;
mproj = zeros(osy,fullwidth_min,'single');
parfor p = 1:ip180
    % load, flat/dark correct, filter, and subpixel shift
    rim0 = zeros(osy,datsize(2),nrings);
    rim180 = rim0;
    for i = 1:nrings
        % load
        im1 = rot90(double(h5read([rawbasedir ringnames{h,i} filesep ringnames{h,i} '.nxs'],h5ImagePath,...
        [1,ycrop(1),p],[datsize(2),osy,1])))-mddarks_crop(:,:,i);
        % flat correct
        im1 = im1./mflats_crop(:,:,i);
        
        % load
        im2 = rot90(double(h5read([rawbasedir ringnames{h,i} filesep ringnames{h,i} '.nxs'],h5ImagePath,...
        [1,ycrop(1),p+ip180],[datsize(2),osy,1])))-mddarks_crop(:,:,i);
        % flat correct
        im2 = im2./mflats_crop(:,:,i);
        
        % zinger correction
        if and(zingerCorrection == 1, not(p==1))
            im1m1 = rot90(double(h5read([rawbasedir ringnames{h,i} filesep ringnames{h,i} '.nxs'],h5ImagePath,...
                [1,ycrop(1),p-1],[datsize(2),osy,1])))-mddarks_crop(:,:,i);
            % flat correct
            im1m1 = im1m1./mflats_crop(:,:,i);
            im2m1 = rot90(double(h5read([rawbasedir ringnames{h,i} filesep ringnames{h,i} '.nxs'],h5ImagePath,...
                [1,ycrop(1),p+ip180-1],[datsize(2),osy,1])))-mddarks_crop(:,:,i);
            % flat correct
            im2m1 = im2m1./mflats_crop(:,:,i);
            
            d1 = im1-im1m1;
            d2 = im2-im2m1;
            
            sig = 3.5;
            fun = @(block_struct) (mean2(block_struct.data)+sig*std2(block_struct.data)) * ones(size(block_struct.data));
            %fun = @(block_struct) prctile(block_struct.data,[99.9],'all') * ones(size(block_struct.data));
            blockSize = [64,64];
            
            B1 = blockproc(d1,blockSize,fun);
            bim1 = d1>B1;
            B2 = blockproc(d2,blockSize,fun);
            bim2 = d2>B1;
            
            im1(bim1) = im1m1(bim1);
            im2(bim2) = im2m1(bim2);
        end
        
        % gain/offset correction
        if gainoffsetCorrection == 1
            im1 = (im1-offsets)./gains;
            im2 = (im2-offsets)./gains;
        end
        
        % subpixel shift
        if i == 1
        else
            im1 = subpixelshift(im1,0,rem(stitchpos(i),1));
            im2 = subpixelshift(im2,0,rem(stitchpos(i),1));
        end
        
        rim0(:,:,i) = im1;
        rim180(:,:,i) = im2;
    end
        
    % % intensity normalizing in overlap regions
    ebump = 10;
    for i = fliplr(1:nrings-1)
        olregi = floor(stitchpos(i+1))+1:floor(stitchpos(i))+datsize(2);
        olregi = datsize(2)-length(olregi)+1:datsize(2);
        olregi = olregi(ebump):olregi(end-ebump);
        olregip1 = (1:length(olregi))+ebump;
        
        mor0i = mean(rim0(:,olregi,i),'all');
        mor0ip1 = mean(rim0(:,olregip1,i+1),'all');
        
        mor180i = mean(rim180(:,olregi,i),'all');
        mor180ip1 = mean(rim180(:,olregip1,i+1),'all');
        
        rim0(:,:,i) = (mor0ip1/mor0i)*rim0(:,:,i);
        rim180(:,:,i) = (mor180ip1/mor180i)*rim180(:,:,i);
    end
    proj = zeros(osy,halfwidth,nrings);
    projp180 = zeros(osy,halfwidth,nrings);
    for i = 1:nrings
        proj(:,floor(stitchpos(i))+1:floor(stitchpos(i))+datsize(2),i) = rim0(:,:,i);
        projp180(:,floor(stitchpos(i))+1:floor(stitchpos(i))+datsize(2),i) = rim180(:,:,i);
    end
    halfproj = BlendProjection(proj,halfmask); % note 10.03.2021: BlendProjection should be improved
    halfprojp180 = BlendProjection(projp180,halfmask);
    halfproj = fliplr(halfproj);
    morl = mean(halfproj(:,halfwidth-ceil(corpix_subpix)+1:end),'all');
    morr = mean(halfprojp180(:,1:ceil(corpix_subpix)),'all');
    halfproj = (morr/morl)*halfproj;
    
    tmpproj = zeros(osy,fullwidth,2);
    tmpproj(:,1:halfwidth,1) = halfproj;
    tmpproj(:,halfwidth-ceil(corpix_subpix)+1:end,2) = halfprojp180;
    fullproj = (BlendProjection(tmpproj,fullmask));
    
    % 
    bump = 5;
    fullproj = fullproj(:,bump+1:end-bump); % knock off edges
    fullproj = padarray(fullproj,[0,bump],'symmetric','both');
    
    % crop to uniform width
    if crop_left > crop_right
        fullproj = subpixelshift(fullproj, 0, 0.5);
    end
    fullproj = fullproj(:,crop_left+1:fullwidth-crop_right);

    % add to mean projection
    mproj = mproj + (fullproj/ip180);
    
    % write to file
    fullproj = single(fullproj);
    if imagefmt == "tiff"
        t = Tiff([writedir 'proj_uf_h' num2str(th) '_p' num2str(p,'%04d') '.tif'], 'w');
        t.setTag(tagstruct); t.write(fullproj); t.close();
    else
        fullproj_filename = [writedir 'proj_uf_h' num2str(th) '_p' ...
            num2str(p,'%04d') '.h5'];
        h5create(fullproj_filename, '/proj', size(fullproj), ...
            Datatype='single');
        h5write(fullproj_filename, '/proj', fullproj);
    end
end
toc
mproj_m150 = medfilt2(mproj,[150,150]);
rproj = -(mproj-mproj_m150);
rproj = single(rproj); % ring profile
mproj = single(mproj); % mean projection
if imagefmt == "tiff"
    t = Tiff([writedir 'rproj_h' num2str(th) '.tif'], 'w');
    t.setTag(tagstruct); t.write(rproj); t.close();
    t = Tiff([writedir 'mproj_h' num2str(th) '.tif'], 'w');
    t.setTag(tagstruct); t.write(mproj); t.close();
else
    rproj_filename = [writedir 'rproj_h' num2str(th) '.h5'];
    h5create(rproj_filename, '/proj', size(rproj), Datatype='single');
    h5write(rproj_filename, '/proj', rproj);
    mproj_filename = [writedir 'mproj_h' num2str(th) '.h5'];
    h5create(mproj_filename, '/proj', size(mproj), Datatype='single');
    h5write(mproj_filename, '/proj', mproj);
end
end
%%
end




