%% Mosaic reconstruction for "mouse6_perf_eth"

%% Toolboxes
addpath(genpath('./utils'))

%% Useful functions
cropaboutcenter = @(im,roisize) im(floor(size(im,1)/2)-floor(roisize(2)/2)+1:floor(size(im,1)/2)+floor(roisize(2)/2),...
    floor(size(im,2)/2-roisize(1)/2)+1:floor(size(im,2)/2+roisize(1)/2));
cropto = @(im,roi) im(roi(3):roi(4),roi(1):roi(2));
cropto3d = @(im,roi) im(roi(3):roi(4),roi(1):roi(2),roi(5):roi(6));
croptostack = @(im,roi) im(roi(3):roi(4),roi(1):roi(2),:);
prepimage = @(im,vr) uint8(255*(double(im)-vr(1))/(vr(2)-vr(1)));
prepimagei16 = @(im,vr) uint16((2^16-1)*(double(im)-vr(1))/(vr(2)-vr(1)));

%% Param file
paramfile = './example/param_files/mouse6_perf_eth_full.txt';

%% Read a few useful variables from the param file
fid = fopen(paramfile);
infoStruct = textscan(fid, '%s %s','Delimiter','\t','CommentStyle','//');
infoStruct = cell2struct(infoStruct{2},infoStruct{1},size(infoStruct,1));
samplename = infoStruct.samplename;
projdir = infoStruct.projpath;
tmp = split(infoStruct.heightstep,'-');
hs_range = str2double(tmp{1}):str2double(tmp{2});
nhs = length(hs_range);
infofile = infoStruct.infofile;
rawdatapath = infoStruct.rawdatapath;
h5ImagePath = infoStruct.h5ImagePath;
h5AnglePath = infoStruct.h5AnglePath;
filtertag = infoStruct.projfiltertype;
filterwidth = str2double(infoStruct.projfilterwidth);
if strcmp(filtertag,'paganin')
    photonEnergy = str2double(infoStruct.photonenergy); % entered into param file as keV
    det_dist_mm = str2double(infoStruct.detdist);
    pixsize_mm = ReadPixelSize_ParamFile(paramfile);
    lambda = lambda_from_E(photonEnergy*1e3); % [m]
    lambda_mm = lambda*1e3;
    deltabeta = filterwidth;
end
if strcmp(filtertag,'paganin')
    filtfunc = @(proj) paganin_filter_stack(proj,filterwidth,pixsize_mm,lambda_mm,det_dist_mm);
elseif strcmp(filtertag,'gauss')
    filtfunc = @(proj) imgaussfilt(proj,filterwidth);
else
    filtfunc = @(proj) proj;
end
fclose(fid);
[nrings,ringnames] = ReadRingsAndNames(infofile,hs_range);

pixsize_mm = ReadPixelSize_ParamFile(paramfile);
[angles,ip180] = ReadAngles_ParamFile(paramfile,h5AnglePath);

%% Path to python scripts
pythonscript_fullpath = '/media/bmc/_home3/handover/MosaicReconstruction/utils/SingleGridrecReconstruction.py'; % full path to your python script "SingleGridrecReconstruction.py"

%% Reconstruct a test slice

%  manually select values to be used for all height steps
manstitchpos = [
        207.6    1848.3    1844.6    1848.3;  % 1
        205.9    1848.3    1848.3    1848.3;  % 2
        205.3    1848.6    1848.6    1847.7;  % 3
        205.4    1848.1    1848.6    1847.8;  % 4
        205.0    1848.3    1848.3    1848.3;  % 5
        205.0    1848.5    1846.0    1847.6;  % 6
        203.6    1848.3    1848.3    1848.3;  % 7
        204.6    1848.3    1848.3    1848.3;  % 8
    ];


testdir = [projdir samplename filesep 'test_slices' filesep];
if not(isfolder(testdir)); mkdir(testdir); end

roi = [14900,14900];

% select heightstep, process projections around y values, reconstruct
h = 3;
ylist = 400:200:1800;

for y = ylist
    ycrop = y-7:y+8;
    % build projections using above stitch positions for given height step
    projsavedir = ProjectionProcessingManualEntry(paramfile,h,...
        manstitchpos(h,:),...
        ycrop,filterwidth,filtertag,'pixsize_mm',pixsize_mm,...
        'lambda_mm',lambda_mm,'det_dist_mm',det_dist_mm);
    [sinoblock] = LoadCroppedSinoblock(projsavedir);
    % ring correction:
    mproj = mean(sinoblock,3);
    rproj = single(mproj-imgaussfilt(mproj,50,'Padding','symmetric')); 
    rproj(abs(rproj)>0.1) = 0;
    % reconstruct
    sinoblock = sinoblock-rproj;
    fullsino = squeeze(sinoblock(round(end/2),:,:));
    fullsino = fullsino(5:end-5,:);
    reco = SingleGridrecReconstruction([testdir ...
        'reco_hs' num2str(h) '_' num2str(y, "%04d")],...
        fullsino,angles,pixsize_mm,...
        pythonscript_fullpath);
end