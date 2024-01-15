function [writedir] = ProjectionProcessingManualOverlap(paramfile,hs,ycrop)

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
if not(isfolder(procdir)); mkdir(procdir); end
if not(isfolder([procdir 'proj_manualoverlap' filesep])); mkdir([procdir 'proj_manualoverlap' filesep]); end
%% 0.5 Load measurement info
% read rings and associated scan names
[nrings,ringnames] = ReadRingsAndNames(infofile,hs);

% read sizes, angles, and .par file info
datsize = ReadScanSize([rawbasedir ringnames{1} filesep ringnames{1} '.nxs'],h5ImagePath);

[angles,ip180] = ReadAngles([rawbasedir ringnames{1} filesep ringnames{1} '.nxs'],h5AnglePath);
ip360 = 2*ip180;

pixsize_mm = ReadPixelSize([rawbasedir ringnames{1} filesep ringnames{1} '.par']);

osy = length(ycrop); % output size y
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
writedir = [procdir 'proj_manualoverlap' filesep];
if not(isfolder(writedir)); mkdir(writedir); end

% tiff settings ()
tagstruct.ImageLength = osy; % y
tagstruct.ImageWidth = datsize(1); % x
tagstruct.BitsPerSample = 32; % single precission
tagstruct.RowsPerStrip = stripheight; % strip size (faster loading of strips)
tagstruct.SamplesPerPixel = 1;
tagstruct.Compression = Tiff.Compression.None;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

% pre-crop flats/darks to avoid repeating crop ip180 times
mddarks_crop = mddarks(ycrop,:,:);
mflats_crop = mflats(ycrop,:,:)-mddarks_crop; % avoid repeating subtraction

% start filtering, blending, stitching
for i = 1:nrings
    fprintf('Stitching and saving projections...\n'); 
    fprintf(['Ring ' num2str(i) '/' num2str(nrings) '\n']); 
    tic;
    mproj = zeros(osy,2048,'single');
    parfor p = 1:ip360
        % load
        im1 = rot90(double(h5read([rawbasedir ringnames{i} filesep ringnames{i} '.nxs'],h5ImagePath,...
        [1,ycrop(1),p],[datsize(2),osy,1])))-mddarks_crop(:,:,i);
        % flat correct
        im1 = single(im1./mflats_crop(:,:,i));
        
        t = Tiff([writedir 'proj_r' num2str(i) '_' num2str(p,'%04d') '.tif'], 'w');
        t.setTag(tagstruct); t.write(im1); t.close();
        
        mproj = mproj + (im1/ip360);
    end
    mproj = single(mproj); % mean projection
    t = Tiff([writedir 'mrean_proj_r' num2str(i) '.tif'], 'w');
    t.setTag(tagstruct); t.write(mproj); t.close();
    toc
end

end

