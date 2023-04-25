function [writedir] = ProjectionProcessingManualEntry(paramfile,hs,manstitchpos,ycrop,filterwidth,filtertag,varargin)
%% Input: Dataset to process
%paramfile = '/media/griffin/_home1/anatomix_feb2021/stitching_reco/examples/braini_paraffin_reembed_hs8_param.txt';
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
%% 0.3 Sort out filtering options
if strcmp(filtertag,'paganin')
    def_pixsize_mm = 6.5e-4; % [mm]
    def_lambda_mm = lambda_from_E(27e3); % [mm]
    def_det_dist_mm = 50; % [mm]
    fprintf('Paganin filtering with delta/beta = %d \n',filterwidth)
    if isempty(varargin)
        pixsize_mm = def_pixsize_mm;
        fprintf('NOTE: Default pixel size of %d mm will be used \n',pixsize_mm)
        lambda_mm = def_lambda_mm;
        fprintf('NOTE: Default wavelength of %d mm will be used \n',lambda_mm)
        det_dist_mm = def_det_dist_mm;
        fprintf('NOTE: Default propagation distance of %d mm will be used \n',det_dist_mm)
    else
    ptag = 0; ltag = 0; dtag = 0;
    for i = 1:length(varargin)
        if strcmp(varargin{i},'pixsize_mm')
            pixsize_mm = varargin{i+1};
            ptag = 1;
        elseif strcmp(varargin{i},'lambda_mm')
            lambda_mm = varargin{i+1};
            ltag = 1;
        elseif strcmp(varargin{i},'det_dist_mm')
            det_dist_mm = varargin{i+1};
            dtag = 1;
        end
    end
    
    if ptag == 1
        fprintf('Using pixel size of %d mm \n',pixsize_mm)
    else
        pixsize_mm = def_pixsize_mm;
        fprintf('NOTE: Default pixel size of %d mm will be used \n',pixsize_mm)
    end
    if ltag == 1
        fprintf('Using wavelength of %d mm \n',lambda_mm)
    else
        lambda_mm = def_lambda_mm;
        fprintf('NOTE: Default wavelength of %d mm will be used \n',lambda_mm)
    end
    if dtag == 1
        fprintf('Using propagation distance of %d mm \n',det_dist_mm)
    else
        det_dist_mm = def_det_dist_mm;
        fprintf('NOTE: Default propagation distance of %d mm will be used \n',det_dist_mm)
    end
    end
elseif strcmp(filtertag,'gauss')
    fprintf('Gaussian filtering with sigma = %d pixels \n', filterwidth)    
else
    fprintf('No filtering will be applied, if you want filtering use ''gauss'' or ''paganin'' \n')
end

if strcmp(filtertag,'paganin')
    filtfunc = @(proj) paganin_filter_stack(proj,filterwidth,pixsize_mm,lambda_mm,det_dist_mm);
elseif strcmp(filtertag,'gauss')
    filtfunc = @(proj) imgaussfilt(proj,filterwidth);
else
    filtfunc = @(proj) proj;
end
%% 0.4 Set up directories
procdir = [procbasedir samplename filesep 'tests_hs' num2str(hs,'%02d') filesep];

if not(isfolder(procdir)); mkdir(procdir); end
if not(isfolder([procdir 'stitched_proj_filtered' filesep])); mkdir([procdir 'stitched_proj_filtered' filesep]); end
%% 0.5 Load measurement info
% read available info from table
[nrings,ringnames] = ReadRingsAndNames(infofile,hs);

% read sizes, angles, and .par file info
datsize = ReadScanSize([rawbasedir ringnames{1} filesep ringnames{1} '.nxs'],h5ImagePath);

[angles,ip180] = ReadAngles([rawbasedir ringnames{1} filesep ringnames{1} '.nxs'],h5AnglePath);
% % commented out by Mattia: it's deceptive if the ycrop passed to the
% % function is not respected
% % also this is a logical fallacy: the ycrop from the param file is
% % ignored in both cases
% if not(isfield(infoStruct,'ycrop'))
% ycrop = 1:datsize(1); 
% end
osy = length(ycrop); % output size y
%% 0.6 Load overlap positions
%tmpdir = [procdir 'parameters' filesep];
%load([tmpdir 'stitching_parameters.mat'],'ol_final_subpix','corpix_subpix','fullwidth','stitchpos','halfwidth')
corpix_subpix = manstitchpos(1);
ol_final_subpix = manstitchpos(2:end);
stitchpos = [0,cumsum(ol_final_subpix)];
halfwidth = ceil(datsize(2)+sum(ol_final_subpix));
fullwidth = halfwidth*2-ceil(corpix_subpix);
%% 1.1 Flats/darks loading
% using mean flat and median dark
mflats = zeros(datsize(1),datsize(2),nrings);
mddarks = zeros(datsize(1),datsize(2),nrings);
for i = 1:nrings
    tmp = rot90(h5read([rawbasedir ringnames{i} filesep 'pre_ref.nxs'],'/ref_0'));
    tmp1 = rot90(h5read([rawbasedir ringnames{i} filesep 'post_ref.nxs'],'/ref_0'));
    flats = double(cat(3,tmp,tmp1)); clear tmp tmp1
    mflats(:,:,i) = mean(flats,3);
    
    darks = rot90(h5read([rawbasedir ringnames{i} filesep 'post_dark.nxs'],'/dark_0'));
    mddarks(:,:,i) = double(median(darks,3));
end
clear darks flats tmp1 tmp
%% 1.2 Flat/dark correction, filtering, stitching, saving
% directory for writing projections:
writedir = [procdir 'stitched_proj_filtered' filesep];

% tiff settings ()
tagstruct.ImageLength = osy; % y
tagstruct.ImageWidth = fullwidth; % x
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

% start filtering, blending, stitching
fprintf('Stitching and saving projections...\n'); tic;
mproj = zeros(osy,fullwidth,'single');
parfor p = 1:ip180
    % load, flat/dark correct, filter, and subpixel shift
    rim0 = zeros(osy,datsize(2),nrings);
    rim180 = rim0;
    for i = 1:nrings
        % load
        im = rot90(double(h5read([rawbasedir ringnames{i} filesep ringnames{i} '.nxs'],h5ImagePath,...
            [1,ycrop(1),p],[datsize(2),osy,1])))-mddarks_crop(:,:,i);
        % flat correct
        im = im./mflats_crop(:,:,i);
        % subpixel shift
        if i == 1
        else
            im = subpixelshift(im,0,rem(stitchpos(i),1));
        end
        rim0(:,:,i) = im;
        
        % load
        im = rot90(double(h5read([rawbasedir ringnames{i} filesep ringnames{i} '.nxs'],h5ImagePath,...
            [1,ycrop(1),p+ip180],[datsize(2),osy,1])))-mddarks_crop(:,:,i);
        % flat correct
        im = im./mflats_crop(:,:,i);
        % subpixel shift
        if i == 1
        else
            im = subpixelshift(im,0,rem(stitchpos(i),1));
        end
        rim180(:,:,i) = im;
    end
    
    % % intensity matching in overlap regions
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
    halfproj = BlendProjection(proj,halfmask);
    halfprojp180 = BlendProjection(projp180,halfmask);
    halfproj = fliplr(halfproj);
    morl = mean(halfproj(:,halfwidth-ceil(corpix_subpix)+1:end),'all');
    morr = mean(halfprojp180(:,1:ceil(corpix_subpix)),'all');
    halfproj = (morr/morl)*halfproj;
    
    tmpproj = zeros(osy,fullwidth,2);
    tmpproj(:,1:halfwidth,1) = halfproj;
    tmpproj(:,halfwidth-ceil(corpix_subpix)+1:end,2) = halfprojp180;
    fullproj = (BlendProjection(tmpproj,fullmask));
    
    % filter stitched projection
    bump = 5;
    fullproj = fullproj(:,bump+1:end-bump); % knock off edges
    fullproj = filtfunc(fullproj);
    fullproj = padarray(fullproj,[0,bump],'symmetric','both');
    
    % add to mean projection
    mproj = mproj + (fullproj/ip180);
    
    % write to tiff file
    fullproj = single(fullproj);
    t = Tiff([writedir 'proj_filtered_' num2str(p,'%04d') '.tif'], 'w');
    t.setTag(tagstruct); t.write(fullproj); t.close();
end
toc
mproj_m150 = medfilt2(mproj,[150,150]);
rproj = -(mproj-mproj_m150);
rproj = single(rproj); % ring profile
mproj = single(mproj); % mean projection
t = Tiff([writedir 'ring_proj.tif'], 'w');
t.setTag(tagstruct); t.write(rproj); t.close();
t = Tiff([writedir 'mean_filtered_proj.tif'], 'w');
t.setTag(tagstruct); t.write(mproj); t.close();
%%
end




