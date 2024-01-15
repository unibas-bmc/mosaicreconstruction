function [writedir] = ProjectionProcessingManualOverlap(paramfile,ycrop)

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
if not(isfolder([procdir 'proj_manualoverlap' filesep])); mkdir([procdir 'proj_manualoverlap' filesep]); end
%% 0.5 Load measurement info
% read available info from table
T = readtable(infofile); c = 0;
pixsize_um = T.pixsize_um(1);
pixsize = pixsize_um*1e-6; % [m]
pixsize_mm = pixsize_um*1e-3;
for i = 1:size(T,1)
    if T.heightstep(i) == hs
        c = c+1;
        scandirs{c} = T.scandir{i};
        disknames{c} = T.diskname{i};
        this_suffix = T.suffix(i);
        if this_suffix == 0
            scansuffixes{c} = '';
        else
            scansuffixes{c} = ['_' num2str(this_suffix)];
        end
        if isnumeric(T.addstring)
            addstring{c} = '';
        else
            addstring{c} = T.addstring{i};
        end
        try
            replaceFlats(c) = T.replaceFlats(i);
            replaceDarks(c) = T.replaceDark(i);
        catch
            replaceFlats(c) = NaN;
            replaceDarks(c) = NaN;
        end
    end
end
nrings = c;
for i = 1:nrings
    ringnames{i} = [scandirs{i} scansuffixes{i} addstring{i}];
end

% read sizes, angles, and .par file info
tmp = h5info([rawbasedir ringnames{1} filesep ringnames{1} '.nxs'],'/flyscan_0001/scan_data/orca_image');
datsize = tmp.Dataspace.Size; % size of images (usually 2048x2048)
angles = abs(h5read([rawbasedir ringnames{1} filesep ringnames{1} '.nxs'],'/flyscan_0001/scan_data/tomo1Rz1'));
ip360 = findindex(angles,360+angles(1));% index of projection with angle 360 degrees
ip180 = floor(ip360/2); % index of projection with angle 180 degrees
angles = angles(1:ip180);
if not(isfield(infoStruct,'ycrop'))
ycrop = 1:datsize(1); 
end
osy = length(ycrop); % output size y
%% 1.1 Flats/darks loading
% using mean flat and median dark
mflats = zeros(datsize(1),datsize(2),nrings);
mddarks = zeros(datsize(1),datsize(2),nrings);
for i = 1:nrings
    try
        tmp = rot90(h5read([rawbasedir ringnames{i} filesep 'pre_ref.nxs'],'/ref_0'));
        tmp1 = rot90(h5read([rawbasedir ringnames{i} filesep 'post_ref.nxs'],'/ref_0'));
        flats = double(cat(3,tmp,tmp1)); clear tmp tmp1
        mflats(:,:,i) = mean(flats,3);
    catch
        % flats are not ok
        if ~isnan(replaceFlats(i))
            % read replacement as specified in excel sheet
            if replaceFlats(i)>0
                replaceScanSuffix = ['_' num2str(replaceFlats(i))];
            else
                replaceScanSuffix = '';
            end
            % same addstring?
            replaceByRingname = [scandirs{i} replaceScanSuffix addstring{i}];
            tmp = rot90(h5read([rawbasedir replaceByRingname filesep 'pre_ref.nxs'],'/ref_0'));
            tmp1 = rot90(h5read([rawbasedir replaceByRingname filesep 'post_ref.nxs'],'/ref_0'));
            flats = double(cat(3,tmp,tmp1)); clear tmp tmp1
            mflats(:,:,i) = mean(flats,3);
        else
            warning(['Cannot read flat field files ' rawbasedir ringnames{i} filesep '*_ref.nxs ' newline 'Please add different suffix in replaceFlats column of excel sheet'])
            keyboard
        end
    end
    try
        darks = rot90(h5read([rawbasedir ringnames{i} filesep 'post_dark.nxs'],'/dark_0'));
        mddarks(:,:,i) = double(median(darks,3));
    catch
        % darks are not ok
        if ~isnan(replaceDarks(i))
            % read replacement as specified in excel sheet
            if replaceDarks(i)>0
                replaceScanSuffix = ['_' num2str(replaceDarks(i))];
            else
                replaceScanSuffix = '';
            end
            % same addstring?
            replaceByRingname = [scandirs{i} replaceScanSuffix addstring{i}];
            darks = rot90(h5read([rawbasedir replaceByRingname filesep 'post_dark.nxs'],'/dark_0'));
            mddarks(:,:,i) = double(median(darks,3));
        else
            warning(['Cannot read dark field file ' rawbasedir ringnames{i} filesep 'post_dark.nxs ' newline 'Please add different suffix in replaceDark column of excel sheet'])
            keyboard
        end
    end
end
clear darks flats tmp1 tmp
%% 1.2 Flat/dark correction, filtering, stitching, saving
% directory for writing projections:
writedir = [projpath samplename '_hs' num2str(hs,'%02d') filesep 'proj_manualoverlap' filesep];
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
        im1 = rot90(double(h5read([rawbasedir ringnames{i} filesep ringnames{i} '.nxs'],'/flyscan_0001/scan_data/orca_image',...
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

