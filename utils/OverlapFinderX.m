function [writedir] = OverlapFinderX(paramfile)
%% Overlap finder
%% Input: Dataset to process
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
%filtfunc = @(im) im; % no filter

% less imporrtant
cor_srange1 = 125; % coarse search range
cor_srange2 = 10; % coarse search range
cor_srange2_step = 0.1;
projtosearch = 10; % # projections to look for stitch position in
olsearchrangecoarse = 50; % search range for stitch pos about motor pos
olsearchrangefine = 10;
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
for h = 1:nhs
th = hs_range(h);
fprintf('Working on height %d (%d/%d)\n',th,h,nhs)
thisRingnames = ringnames(h,:);
txpos = txpos_all(h,:);
%% 1.1 Flats/darks loading
% using mean flat and median dark
mflats = zeros(datsize(1),datsize(2),nrings);
mddarks = zeros(datsize(1),datsize(2),nrings);
for i = 1:nrings
    tmp = rot90(h5read([rawbasedir thisRingnames{i} filesep 'pre_ref.nxs'],'/ref_0'));
    tmp1 = rot90(h5read([rawbasedir thisRingnames{i} filesep 'post_ref.nxs'],'/ref_0'));
    flats = double(cat(3,tmp,tmp1)); clear tmp tmp1
    mflats(:,:,i) = mean(flats,3);
    darks = rot90(h5read([rawbasedir thisRingnames{i} filesep 'post_dark.nxs'],'/dark_0'));
    mddarks(:,:,i) = double(median(darks,3));
end

clear darks flats tmp1 tmp
toc
%% 1.2 Find COR and stitch positions
fprintf(['Loading projections at 0 and 180...\n']); tic
% % COR:
cor_guess = 2*round((datsize(2)/2)+(txpos(1)/pixsize_mm)); % guess COR from motor position

% find a projection with content in the window
ccwin = 1:cor_guess; % window to compare 0 and 180
prange = floor(linspace(1,ip180,projtosearch));
contval = zeros(1,length(prange));
for p = 1:length(prange)
    tmp = rot90(double(h5read([rawbasedir thisRingnames{1} filesep thisRingnames{1} '.nxs'],h5ImagePath,...
        [ccwin(1),1,prange(p)],[ccwin(end),datsize(2),1])))-mddarks(:,ccwin,1); % load image at 0 deg
    im0 = tmp./(mflats(:,ccwin,1)-mddarks(:,ccwin,1)); % flat corrrected 0 deg image
    contval(p) = entropy(im0);
end

% load 0 and 180
p = prange(findindex(contval,max(contval)));
tmp = rot90(double(h5read([rawbasedir thisRingnames{1} filesep thisRingnames{1} '.nxs'],h5ImagePath,...
        [1,1,p],[datsize(1),datsize(2),1])))-mddarks(:,:,1); % load image at 0 deg
im0 = filtfunc(tmp./(mflats(:,:,1)-mddarks(:,:,1))); % flat corrrected 0 deg image
tmp = rot90(double(h5read([rawbasedir thisRingnames{1} filesep thisRingnames{1} '.nxs'],h5ImagePath,...
    [1,1,p+ip180],[datsize(1),datsize(2),1])))-mddarks(:,:,1); % load image at 180 deg
tmp = tmp./(mflats(:,:,1)-mddarks(:,:,1)); % flat corrected image at 180
im180 = filtfunc(fliplr(tmp)); % mirror image at 180
if verboseMode
    if batchStartupOptionUsed
        fh = figure('visible','off');
    else
        fh = figure();
    end
    subplot(121), imagesc(im180,prctile(im180,[2,98],'all')'), axis equal tight, title('proj at 180')
    subplot(122), imagesc(im0,prctile(im180,[2,98],'all')'), axis equal tight, title('proj at 0')
    print(fh,'-dpng',[figdir 'h' num2str(h) '_images_0_180deg.png'])
end
toc
% coarse search
ccwin = 1:cor_guess/2;%cor_guess; % window to compare 0 and 180
scanRange = max(cor_guess-cor_srange1,1):(cor_guess+cor_srange1); % number of overlap pixels. Range over which to search for best overlap
ccval = zeros(1,length(scanRange));
fprintf(['Finding center of rotation (coarse)...\n']); tic
for il = 1:length(scanRange)
    tim180 = subpixelshift(im180,0,scanRange(il));
    ccval(il) = corr2(im0(:,ccwin),...
        tim180(:,ccwin));
end; toc
corpix = scanRange(ccval == max(ccval));
if verboseMode
    if batchStartupOptionUsed
        fh = figure('visible','off');
    else
        fh = figure();
    end
    subplot(121), plot(scanRange,ccval,'.-'), vline(corpix,'r:',num2str(corpix))
    vline(cor_guess,'b:','motor pos.')
    xlabel('x-shift [pix]'), ylabel('cc(T_x[180],0)'), axis tight
end

% refined search
ccwin = 1:corpix; % window to compare 0 and 180
scanRange = max(corpix-cor_srange2,1):cor_srange2_step:corpix+cor_srange2; % search about previous max
ccval = zeros(1,length(scanRange));
fprintf(['Finding center of rotation (refining)...\n']); tic
for il = 1:length(scanRange)
    tim180 = subpixelshift(im180,0,scanRange(il));
    ccval(il) = corr2(im0(:,ccwin),...
        tim180(:,ccwin));
end; toc
corpix = round(scanRange(ccval == max(ccval)));
corpix_subpix = scanRange(ccval == max(ccval));
if verboseMode
    subplot(122), plot(scanRange,ccval,'.-')
    vline(cor_guess,'b:','motor pos.')
    vline(corpix_subpix,'k:',num2str(corpix_subpix))
    vline(corpix,'r:',num2str(corpix))
    xlabel('x-shift [pix]'), ylabel('cc(T_x[180],0)'), axis tight
    print(fh,'-dpng',[figdir 'h' num2str(h) '_centerOfRotation.png'])
end

% % Stitching positions
fprintf(['Finding ring stitching positions...\n']);
% get stitch position estimates from the motors
pixshifts = zeros(1,nrings-1);
for i = 1:nrings-1
    this_xshift = txpos(i)-txpos(i+1);
    pixshifts(i) = round(this_xshift/pixsize_mm);
end

% search a range of projections for stitch positions
prange = linspace(findindex(angles,anglestosearch(1)),findindex(angles,anglestosearch(2)),projtosearch); % subset of projections to find stitching positions
prange = round(prange); % range of projections to look in
ccvals = zeros(length((-ceil(olsearchrangefine/2)+1):olsearchstep:floor(olsearchrangefine/2)),length(prange),nrings-1);
scanRanges = ccvals;
ols = zeros(length(prange),nrings-1);
promsmat = zeros(length(prange),nrings-1);
locsmat = zeros(length(prange),nrings-1);
promsmat1 = zeros(length(prange),nrings-1);
locsmat1 = zeros(length(prange),nrings-1);
if verboseMode
    if batchStartupOptionUsed
        fh = figure('visible','off');
    else
        fh = figure();
    end
end
for i = 1:nrings-1
    fprintf(['Finding overlap of rings ' num2str(i) ' and ' num2str(i+1) '...\n'])
    tic
    for p = 1:length(prange)
    
    tmp = rot90(double(h5read([rawbasedir thisRingnames{i} filesep thisRingnames{i} '.nxs'],h5ImagePath,...
        [1,1,prange(p)],[datsize(1),datsize(2),1])))-mddarks(:,:,i);
    imr1 = filtfunc(tmp./(mflats(:,:,i)-mddarks(:,:,i)));
    tmp = rot90(double(h5read([rawbasedir thisRingnames{i+1} filesep thisRingnames{i+1} '.nxs'],h5ImagePath,...
        [1,1,prange(p)],[datsize(1),datsize(2),1])))-mddarks(:,:,i+1);
    imr2 = filtfunc(tmp./(mflats(:,:,i+1)-mddarks(:,:,i+1)));
    
    ccwin1 = datsize(2)-windowwidth+1:datsize(2);    
    imr1c = imr1(:,ccwin1);
    
    quickRange = pixshifts(i)+((-ceil(olsearchrangecoarse/2)+1):floor(olsearchrangecoarse/2));
    quickccval = zeros(1,length(quickRange));
    for il = 1:length(quickRange)
        timr2c = subpixelshift(imr2,0,quickRange(il));
        timr2c = timr2c(:,ccwin1);
        %timr2c = imr2(:,ccwin1-scanRange(il));
        quickccval(il) = corr2(imr1c,timr2c);
    end
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
    for il = 1:length(scanRange)
        timr2c = subpixelshift(imr2,0,scanRange(il));
        timr2c = timr2c(:,ccwin1);
        ccval(il) = corr2(imr1c,timr2c);
    end
    this_ol = scanRange(ccval==max(ccval));
    this_ol = this_ol(1);
    scanRanges(:,p,i) = scanRange;
    ccvals(:,p,i) = ccval;
    ols(p,i) = this_ol;
    
    [~,locs,~,proms] = findpeaks(ccval);
    [~,I] = sort(proms,'descend');
    if ne(length(I),0)
    locsvec = locs(I(1));
    promsvec = proms(I(1));
    else
    locsvec = round(length(scanRange)/2); 
    promsvec = 0;
    end
    promsmat(p,i) = promsvec;
    locsmat(p,i) = locsvec;
    
    end
    if verboseMode
        subplot(nrings-1,1,i)
        plot(scanRanges(:,:,i),ccvals(:,:,i),'.-'), hold on
        plot(ols(:,i),max(ccvals(:,:,i)),'r*')
        vline(mean(ols(:,i),1),'r-.',num2str(mean(ols(:,i),1)))
        title(['rings ' num2str(i) ' and ' num2str(i+1)]), xlabel('x-shift [pix]'), ylabel('cc(T_x[180],0)')
        if ~batchStartupOptionUsed
            drawnow
        end
    end
    toc
end
if verboseMode
    print(fh,'-dpng',[figdir 'h' num2str(h) '_ringOverlaps.png'])
end

if verboseMode
for r = 1:nrings-1
fh = figure(); subplot(211)
plot(promsmat(:,r),'.-'), xlabel('proj #'), ylabel('peak prominence')
subplot(212)
plot(locsmat(:,r),'.-'), xlabel('proj #'), ylabel('peak position')
sgtitle(['rings ' num2str(r) ' and ' num2str(r+1)])
drawnow
print(fh,'-dpng',[figdir 'h' num2str(h) '_ringOverlapPeaks ' num2str(r) '.png'])
end
end

W = promsmat1./sum(promsmat1,1);
ol_final_subpix = sum(W.*ols,1); % prominence-weighted mean
ol_final = round(ol_final_subpix);
stitchpos = [0,cumsum(ol_final_subpix)];
halfwidth = ceil(datsize(2)+sum(ol_final_subpix)); % width before stitching off-axis
fullwidth = halfwidth*2-ceil(corpix_subpix); % width of final stitched projections (off-axis)

%% Output: write information to a file
% % write parameters and found positions to file
writedir = [procdir 'parameters' filesep];
save([writedir 'stitching_parameters_h' num2str(h) '.mat'],'ol_final_subpix','corpix_subpix','fullwidth','stitchpos','halfwidth')

t = table(corpix_subpix,ol_final_subpix,stitchpos,halfwidth,fullwidth);
writetable(t,[writedir 'stitching_parameters_h' num2str(h) '.csv'])
%%
end
