function [writedir] = OverlapFinderY(paramfile)
%% Overlap finder
%% Input: Dataset to process
%paramfile = './example/param_files/mouse6_perf_eth_full.txt';
fprintf(['Initializing, loading flats/darks...\n']); tic
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
projdir = infoStruct.projpath;
verboseMode = logical(str2num(infoStruct.verboseMode)); % shows plots if true

h5ImagePath = infoStruct.h5ImagePath;
h5AnglePath = infoStruct.h5AnglePath;
%% 0.3 Internal settings (for debugging)
% more important
anglestosearch = [0,180]; % angle range to look for stitching in; choose where sample is in full fov
filtfunc = @(im) medfilt2(im,[3,3],'symmetric'); % filter for projection for overlap finder
filtfunc = @(im) imgaussfilt(im,1); % filter for projection for overlap finder
%filtfunc = @(im) im; % no filter

% less imporrtant
projtosearch = 5; % # projections to look for stitch position in
olsearchrangecoarse = 50; % search range for stitch pos about motor pos
olsearchrangefine = 15;
olsearchstep = 0.25;
windowwidth = 250;
%% 0.4 Set up directories
procdir = [projdir samplename filesep];
if not(isfolder(procdir)); mkdir(procdir); end
if not(isfolder([procdir 'parameters' filesep])); mkdir([procdir 'parameters' filesep]); end
if verboseMode 
    figdir = [procdir filesep 'figs' filesep];
    if not(isfolder(figdir)); mkdir(figdir); end
end
readdir = [procdir 'stitched_proj_filtered' filesep];
%% 0.5 Load measurement info
% read rings and associated scan names
[nrings,ringnames] = ReadRingsAndNames(infofile,hs_range);

% read sizes, angles, and .par file info
datsize = ReadScanSize([rawbasedir ringnames{1,1} filesep ringnames{1,1} '.nxs'],h5ImagePath);

[angles,ip180] = ReadAngles([rawbasedir ringnames{1,1} filesep ringnames{1,1} '.nxs'],h5AnglePath);

pixsize_mm = ReadPixelSize([rawbasedir ringnames{1,1} filesep ringnames{1,1} '.par']);

txpos_all = zeros(nhs,nrings); % get the tx positions to guess cor and stitch positions
tzpos_all = zeros(nhs,nrings); % get tz positions to know height step
for io = 1:nhs
for i = 1:nrings
    txpos_all(io,i) = ReadTxPos([rawbasedir ringnames{io,i} filesep ringnames{io,i} '.par']);
    tzpos_all(io,i) = ReadTzPos([rawbasedir ringnames{io,i} filesep ringnames{io,i} '.par']);
end
end
%% 1.0 Loop over heights
fprintf(['Finding height step stitching positions...\n']);
%% 1.1 Find stitch positions
% % Stitching positions
% get stitch position estimates from the motors
tzpos = tzpos_all(:,1);
pixshifts = round(abs(diff(tzpos))/pixsize_mm);

% search a range of projections for stitch positions
prange = linspace(findindex(angles,anglestosearch(1)),findindex(angles,anglestosearch(2)),projtosearch); % subset of projections to find stitching positions
prange = round(prange); % range of projections to look in
quickccvals = zeros(length((-ceil(olsearchrangecoarse/2)+1):floor(olsearchrangecoarse/2)),length(prange),nhs-1);
quickscanRanges = quickccvals;
ccvals = zeros(length((-ceil(olsearchrangefine/2)+1):olsearchstep:floor(olsearchrangefine/2)),length(prange),nhs-1);
scanRanges = ccvals;
ols = zeros(length(prange),nhs-1);
promsmat = zeros(length(prange),nhs-1);
locsmat = zeros(length(prange),nhs-1);
promsmat1 = zeros(length(prange),nhs-1);
locsmat1 = zeros(length(prange),nhs-1);
for hi = 1:nhs-1
h = hs_range(hi);
fprintf('Working on heights %d and %d (%d/%d)\n',h,h+1,hi,nhs-1)
tic
figure
title(['Heights ' num2str(h) ' and ' num2str(h+1)])
xlabel('y overlap [pixel]'), ylabel('cross correlation')
hold on
for p = 1:length(prange)
    imr1 = filtfunc(imread([readdir 'proj_uf_h' num2str(h) '_p' num2str(p,'%04d') '.tif']));
    imr2 = filtfunc(imread([readdir 'proj_uf_h' num2str(h+1) '_p' num2str(p,'%04d') '.tif']));
    imr1 = flipud(imr1);
    imr2 = flipud(imr2);
    
    ccwin1 = datsize(2)-windowwidth+1:datsize(2);
    imr1c = imr1(ccwin1,:);
    
    quickRange = pixshifts(hi)+((-ceil(olsearchrangecoarse/2)+1):floor(olsearchrangecoarse/2));
    quickccval = zeros(1,length(quickRange));
    parfor il = 1:length(quickRange)
        timr2c = subpixelshift(imr2,quickRange(il),0);
        timr2c = timr2c(ccwin1,:);
        quickccval(il) = corr2(imr1c,timr2c);
    end
    quickscanRanges(:,p,hi) = quickRange;
    quickccvals(:,p,hi) = quickccval;
    
    [~,locs,~,proms] = findpeaks(quickccval);
    [~,I] = sort(proms,'descend');
    %tmp_max = quickRange(findindex(quickccval,max(quickccval)));
    if ne(length(I),0)
        locsvec = locs(I(1));
        promsvec = proms(I(1));
    else
        locsvec = round(length(quickRange)/2);
        promsvec = 0;
    end
    promsmat1(p,i) = promsvec;
    locsmat1(p,i) = locsvec;
    locs = locsvec;
    tmp_max = quickRange(locs(1));
    scanRange = tmp_max + ((-ceil(olsearchrangefine/2)+1):olsearchstep:floor(olsearchrangefine/2));
    ccval = zeros(1,length(scanRange));
    parfor il = 1:length(scanRange)
        timr2c = subpixelshift(imr2,scanRange(il),0);
        timr2c = timr2c(ccwin1,:);
        ccval(il) = corr2(imr1c,timr2c);
    end
    this_ol = scanRange(ccval==max(ccval));
    scanRanges(:,p,hi) = scanRange;
    ccvals(:,p,hi) = ccval;
    ols(p,hi) = this_ol;
    
    [~,locs,~,proms] = findpeaks(ccval);
    [~,I] = sort(proms,'descend');
    if ne(length(I),0)
        locsvec = locs(I(1));
        promsvec = proms(I(1));
    else
        locsvec = round(length(scanRange)/2);
        promsvec = 0;
    end
    promsmat(p,hi) = promsvec;
    locsmat(p,hi) = locsvec;
    
    tmp = plot(quickscanRanges(:,p,hi),quickccvals(:,p,hi),'.-','DisplayName',[num2str(prange(p))]);
    plot(scanRanges(:,p,hi),ccvals(:,p,hi),'.-','HandleVisibility','off','Color',tmp.Color)
    drawnow
end
toc
legend
end

W = promsmat./sum(promsmat,1);
manstitchposy = sum(W.*ols,1); % prominence-weighted mean

%% Output: write information to a file
% % write parameters and found positions to file
writedir = [procdir 'parameters' filesep];
save([writedir 'stitching_parameters_y.mat'],'manstitchposy')

t = table(manstitchposy);
writetable(t,[writedir 'stitching_parameters_y.csv'])
%%
end
