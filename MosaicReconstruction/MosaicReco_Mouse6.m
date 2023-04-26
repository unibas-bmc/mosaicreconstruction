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

%% Build a roughly stitched mosaic projection
projNo = 1;
medkernel = [7,7];
gausskernel = 3;
[proj,writedir] = PreviewStitchProjectionXY(paramfile,projNo,medkernel,gausskernel);

%% Automatically find COR and X overlap positions (all height steps)
[writedir] = OverlapFinderX(paramfile);

%% Check selected overlap positions manually
% % load automatically found positions
stitchparamdir = [projdir samplename filesep 'parameters' filesep];
olpix = zeros(nrings-1,nhs);
corpix = zeros(1,nhs);
for h = 1:nhs
    t = readtable([stitchparamdir 'stitching_parameters_h' num2str(h) '.csv']);
    olpix(:,h) = [t.ol_final_subpix_1,t.ol_final_subpix_2,t.ol_final_subpix_3];
    corpix(h) = t.corpix_subpix;
end

figure, plot(olpix(:),'o'); hold on; a = gca;
plot(a.XLim,[1,1]*nanmedian(olpix(:)),'r:','LineWidth',2)
text(a.XLim(1)+diff(a.XLim)/30,nanmedian(olpix(:))-diff(a.YLim)/30,['median: ' num2str(nanmedian(olpix(:)))])
xlabel('scan'), ylabel('estimated stitch position')

figure, plot(corpix,'o'); hold on; a = gca;
plot(a.XLim,[1,1]*nanmedian(corpix),'r:','LineWidth',2)
text(a.XLim(1)+diff(a.XLim)/30,nanmedian(corpix)-diff(a.YLim)/30,['median: ' num2str(nanmedian(corpix))])
xlabel('height step'), ylabel('estimated COR position')

%  manually select values to be used for all height steps
% cor = 205;
% s1x = 1848.3;
% s2x = 1848.3;
% s3x = 1848.3;
% manstitchpos = [cor,s1x,s2x,s3x];
autostitchpos = [corpix', olpix'];
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

% loop over heights, process projections for ycrop, reconstruct
ycrop = 1017:1032;

testdir = [projdir samplename filesep 'stitchpos_tests' filesep];
if not(isfolder(testdir)); mkdir(testdir); end

% note: we do a ring correction here
roi = [14900,14900];
recos = zeros(roi(2),roi(1),nhs,'single');
for h = 1:nhs
    % build projections using above stitch positions for given height step
    projdir = ProjectionProcessingManualEntry(paramfile,h,...
        manstitchpos(h,:),...
        ycrop,filterwidth,filtertag,'pixsize_mm',pixsize_mm,...
        'lambda_mm',lambda_mm,'det_dist_mm',det_dist_mm);
    [sinoblock] = LoadCroppedSinoblock(projdir);
    % ring correction:
    mproj = mean(sinoblock,3);
    rproj = single(mproj-imgaussfilt(mproj,50,'Padding','symmetric')); 
    rproj(abs(rproj)>0.1) = 0;
    % reconstruct
    sinoblock = sinoblock-rproj;
    fullsino = squeeze(sinoblock(round(end/2),:,:));
    fullsino = fullsino(5:end-5,:);
    reco = SingleGridrecReconstruction([testdir ...
        'reco_manualoverlap_hs' num2str(h)],fullsino,angles,pixsize_mm,...
        pythonscript_fullpath);
    recos(:,:,h) = cropaboutcenter(reco,roi);
end

if isfile([testdir 'reco_manualoverlap_allhs.h5']); delete([testdir 'reco_manualoverlap_allhs.h5']); end
h5create([testdir 'reco_manualoverlap_allhs.h5'],'/reco',size(recos),'Datatype','single')
h5write([testdir 'reco_manualoverlap_allhs.h5'],'/reco',recos)

%% Tweak overlap positions manually
% % set up a directory for tests
testdir = [projdir samplename filesep 'stitchpos_tests' filesep];
if not(isfolder(testdir)); mkdir(testdir); end

% % generates a stack of cropped projections before stitching
this_hs = 1;
this_ycrop = 1024-7:1024+8;
readdir = ProjectionProcessingManualOverlap(paramfile,this_hs,this_ycrop);
[projvol,mprojvol] = LoadProjectionsManualOverlap(paramfile,this_hs);

sliceNo = 9;

% % Paganin filter projections
projvol_pag = projvol;
for i1 = 1:size(projvol,4)
    projvol_pag(:,:,:,i1) = filtfunc(projvol(:,:,:,i1));
end

% % check center of rotation
corRange = 205.9-8:205.9+8;
% note: motor position would be cor_guess, found position would be cor_subpix
padSize = 2000;
cropSize = [3000,3000];
recos_cor = zeros(cropSize(2),cropSize(1),length(corRange),'single');
for i1 = 1:length(corRange)
    this_cor = corRange(i1);
    halfwidth = 2048;
    fullwidth = halfwidth*2-ceil(this_cor);
    
    sino1 = flipud(squeeze(projvol_pag(sliceNo,:,1:ip180,1)))';
    sino2 = squeeze(projvol_pag(sliceNo,:,ip180+1:2*ip180,1))';
    sino2 = subpixelshift(sino2,0,this_cor-ceil(this_cor));
    
    
    fullmask = zeros(ip180,fullwidth,2);
    fullmask(:,1:halfwidth-1,1) = ones(ip180,halfwidth-1);
    fullmask(:,halfwidth-ceil(this_cor)+2:end,2) = ones(ip180,halfwidth-1);
    fullsino = zeros(ip180,fullwidth,2);
    fullsino(:,1:2048,1) = sino1;
    fullsino(:,2048-ceil(this_cor)+1:end,2) = sino2;
    
    fullsino = (BlendProjection(fullsino,fullmask))';
    fullsino = padarray(fullsino,[padSize,0],'symmetric','both'); % padding since this is local

    
    reco = SingleGridrecReconstruction([testdir ...
        'localreco_COR_' num2str(this_cor)],fullsino,angles,pixsize_mm,...
        pythonscript_fullpath);
    reco = reco(padSize+1:end-padSize,padSize+1:end-padSize);

    recos_cor(:,:,i1) = cropaboutcenter(reco,cropSize);
end
figure, imshow3D(recos_cor,prctile(recos_cor,[1,99],'all')')

%% check stitching positions of each ring
% % Ring 1
this_cor = 205.9;

olpix1_range = 1850.4-5:1850.4+5;
this_nrings = 2;

padSize = 750;
cropSize = [7500,7500];
recos_crop = zeros(cropSize(2),cropSize(1),length(olpix1_range),'single');
for i1 = 1:length(olpix1_range)
    stitchpos = [0,cumsum(olpix1_range(i1))];
    halfwidth = ceil(2048+sum(stitchpos));
    fullwidth = halfwidth*2-ceil(this_cor);
    
    sino0 = zeros(ip180,2048,this_nrings);
    sino180 = sino0;
    for r = 1:this_nrings
        sino1 = squeeze(projvol_pag(sliceNo,:,1:ip180,r))';
        sino2 = squeeze(projvol_pag(sliceNo,:,ip180+1:2*ip180,r))';
        if r == 1
        else
            sino1 = subpixelshift(sino1,0,rem(stitchpos(r),1));
            sino2 = subpixelshift(sino2,0,rem(stitchpos(r),1));
        end
        sino0(:,:,r) = sino1;
        sino180(:,:,r) = sino2;
    end
    proj = zeros(ip180,halfwidth,this_nrings);
    projp180 = zeros(ip180,halfwidth,this_nrings);
    for r = 1:this_nrings
        proj(:,floor(stitchpos(r))+1:floor(stitchpos(r))+2048,r) = sino0(:,:,r);
        projp180(:,floor(stitchpos(r))+1:floor(stitchpos(r))+2048,r) = sino180(:,:,r);
    end
    halfmask = zeros(ip180,halfwidth,this_nrings);
    for r = 1:this_nrings
        if r == 1
            halfmask(:,floor(stitchpos(r))+1:floor(stitchpos(r))+2048-1,r) = ones(ip180,2048-1);
            tmp = floor(stitchpos(r+1))+1:floor(stitchpos(r))+2048-1;
            halfmask(:,tmp,r) = ones(ip180,length(tmp)).*linspace(1,0,length(tmp));
        elseif r == this_nrings
            halfmask(:,floor(stitchpos(r))+2:floor(stitchpos(r))+2048,r) = ones(ip180,2048-1);
            tmp = floor(stitchpos(r))+1:floor(stitchpos(r-1))+2048-1;
            halfmask(:,tmp,r) = ones(ip180,length(tmp)).*linspace(0,1,length(tmp));
        else
            halfmask(:,floor(stitchpos(r))+2:floor(stitchpos(r))+2048-1,r) = ones(ip180,2048-2);
            tmp = floor(stitchpos(r+1))+1:floor(stitchpos(r))+2048-1;
            halfmask(:,tmp,r) = ones(ip180,length(tmp)).*linspace(1,0,length(tmp));
            tmp = floor(stitchpos(r))+1:floor(stitchpos(r-1))+2048-1;
            halfmask(:,tmp,r) = ones(ip180,length(tmp)).*linspace(0,1,length(tmp));
        end
    end

    halfproj = BlendProjection(proj,halfmask); % note 10.03.2021: blendproj should be improved
    halfprojp180 = BlendProjection(projp180,halfmask);
    halfproj = fliplr(halfproj);
    morl = mean(halfproj(:,halfwidth-ceil(this_cor)+1:end),'all');
    morr = mean(halfprojp180(:,1:ceil(this_cor)),'all');
    halfproj = (morr/morl)*halfproj;
    
    tmpproj = zeros(ip180,fullwidth,2);
    tmpproj(:,1:halfwidth,1) = halfproj;
    tmpproj(:,halfwidth-ceil(this_cor)+1:end,2) = halfprojp180;
    fullmask = zeros(ip180,fullwidth,2);
    fullmask(:,1:halfwidth-1,1) = ones(ip180,halfwidth-1);
    fullmask(:,halfwidth-ceil(this_cor)+2:end,2) = ones(ip180,halfwidth-1);
    tmp = halfwidth-ceil(this_cor)+2:halfwidth;
    fullmask(:,halfwidth-ceil(this_cor)+2:halfwidth,1) = ones(ip180,length(tmp)).*linspace(1,0,length(tmp));
    fullmask(:,halfwidth-ceil(this_cor)+2:halfwidth,2) = ones(ip180,length(tmp)).*linspace(0,1,length(tmp));
    fullproj = (BlendProjection(tmpproj,fullmask))';
    fullproj = fullproj(2:end-1,:);
    
    fullproj = padarray(fullproj,[padSize,0],'symmetric','both'); % padding since this is local

    reco = SingleGridrecReconstruction([testdir ...
        'localreco_ol_r1_' num2str(olpix1_range(i1))],fullproj,angles,pixsize_mm,...
        pythonscript_fullpath);
    reco = reco(padSize+1:end-padSize,padSize+1:end-padSize);

    recos_crop(:,:,i1) = cropaboutcenter(reco,cropSize);
end

%gsrange = prctile(recos_crop,[1,99],'all')';
gsrange = [-0.01,0.075];

figure, imshow3D(recos_crop,gsrange)

cent = [0,0];
rad = 2048-this_cor/2;
ax = -cropSize(1)/2+1:cropSize(1)/2;
figure, imagesc(ax,ax,recos_crop(:,:,round(end/2)),gsrange)
axis equal tight, colormap gray
hold on
%viscircles(cent,rad,'EnhanceVisibility',false);
rectangle('Position',[cent(1)-rad,cent(2)-rad,rad*2,rad*2],'Curvature',[1,1],...
    'EdgeColor','r')

%% Tweak any
cor_range = 205.9;
s1_range = 1848.3-8:1848.3+8;
s2_range = 1848.3;
s3_range = 1848.3;

min_size = ceil(2048+max([0,cumsum([min(s1_range),min(s2_range),min(s2_range)])]))*2-ceil(min(cor_range));

padSize = 0;
cropSize = [14000,14000];
recos_crop = zeros(cropSize(2),cropSize(1),length(cor_range),length(s1_range),...
    length(s2_range),length(s3_range),'single');
for i1 = 1:length(cor_range)
for i2 = 1:length(s1_range)
for i3 = 1:length(s2_range)
for i4 = 1:length(s3_range)
    these_olpix = [s1_range(i2),s2_range(i3),s3_range(i4)];
    this_cor = cor_range(i1);
    stitchpos = [0,cumsum(these_olpix)];
    halfwidth = ceil(2048+max(stitchpos));
    fullwidth = halfwidth*2-ceil(this_cor);
    
    sino0 = zeros(ip180,2048,nrings);
    sino180 = sino0;
    for r = 1:nrings
        sino1 = squeeze(projvol_pag(sliceNo,:,1:ip180,r))';
        sino2 = squeeze(projvol_pag(sliceNo,:,ip180+1:2*ip180,r))';
        if r == 1
        else
            sino1 = subpixelshift(sino1,0,rem(stitchpos(r),1));
            sino2 = subpixelshift(sino2,0,rem(stitchpos(r),1));
        end
        sino0(:,:,r) = sino1;
        sino180(:,:,r) = sino2;
    end
    
    proj = zeros(ip180,halfwidth,nrings);
    projp180 = zeros(ip180,halfwidth,nrings);
    for r = 1:nrings
        proj(:,floor(stitchpos(r))+1:floor(stitchpos(r))+2048,r) = sino0(:,:,r);
        projp180(:,floor(stitchpos(r))+1:floor(stitchpos(r))+2048,r) = sino180(:,:,r);
    end
    halfmask = zeros(ip180,halfwidth,nrings);
    for r = 1:nrings
        if r == 1
            halfmask(:,floor(stitchpos(r))+1:floor(stitchpos(r))+2048-1,r) = ones(ip180,2048-1);
            tmp = floor(stitchpos(r+1))+1:floor(stitchpos(r))+2048-1;
            halfmask(:,tmp,r) = ones(ip180,length(tmp)).*linspace(1,0,length(tmp));
        elseif r == nrings
            halfmask(:,floor(stitchpos(r))+2:floor(stitchpos(r))+2048,r) = ones(ip180,2048-1);
            tmp = floor(stitchpos(r))+1:floor(stitchpos(r-1))+2048-1;
            halfmask(:,tmp,r) = ones(ip180,length(tmp)).*linspace(0,1,length(tmp));
        else
            halfmask(:,floor(stitchpos(r))+2:floor(stitchpos(r))+2048-1,r) = ones(ip180,2048-2);
            tmp = floor(stitchpos(r+1))+1:floor(stitchpos(r))+2048-1;
            halfmask(:,tmp,r) = ones(ip180,length(tmp)).*linspace(1,0,length(tmp));
            tmp = floor(stitchpos(r))+1:floor(stitchpos(r-1))+2048-1;
            halfmask(:,tmp,r) = ones(ip180,length(tmp)).*linspace(0,1,length(tmp));
        end
    end

    halfproj = BlendProjection(proj,halfmask); % note 10.03.2021: blendproj should be improved
    halfprojp180 = BlendProjection(projp180,halfmask);
    halfproj = fliplr(halfproj);
    morl = mean(halfproj(:,halfwidth-ceil(this_cor)+1:end),'all');
    morr = mean(halfprojp180(:,1:ceil(this_cor)),'all');
    halfproj = (morr/morl)*halfproj;
    
    tmpproj = zeros(ip180,fullwidth,2);
    tmpproj(:,1:halfwidth,1) = halfproj;
    tmpproj(:,halfwidth-ceil(this_cor)+1:end,2) = halfprojp180;
    fullmask = zeros(ip180,fullwidth,2);
    fullmask(:,1:halfwidth-1,1) = ones(ip180,halfwidth-1);
    fullmask(:,halfwidth-ceil(this_cor)+2:end,2) = ones(ip180,halfwidth-1);
    tmp = halfwidth-ceil(this_cor)+2:halfwidth;
    fullmask(:,halfwidth-ceil(this_cor)+2:halfwidth,1) = ones(ip180,length(tmp)).*linspace(1,0,length(tmp));
    fullmask(:,halfwidth-ceil(this_cor)+2:halfwidth,2) = ones(ip180,length(tmp)).*linspace(0,1,length(tmp));
    fullproj = (BlendProjection(tmpproj,fullmask))';
    fullproj = fullproj(2:end-1,:);
    
    fullproj = padarray(fullproj,[padSize,0],'symmetric','both'); % padding since this is local

    reco = SingleGridrecReconstruction([testdir ...
        'cor_' num2str(cor_range(i1)) '_s1_' num2str(s1_range(i2)) ...
        '_s2_' num2str(s2_range(i3)) '_s3_' num2str(s3_range(i4)) ],...
        fullproj,angles,pixsize_mm,...
        pythonscript_fullpath);
    reco = reco(padSize+1:end-padSize,padSize+1:end-padSize);

    recos_crop(:,:,i1,i2,i3,i4) = cropaboutcenter(reco,cropSize);
end
end
end
end

gsrange = [-0.01,0.05];

figure, imshow3D(squeeze(recos_crop),gsrange)

cent = [0,0];
rad = [2048-cor_range(1)/2,2048-cor_range(1)/2+s1_range(1),...
    2048-cor_range(1)/2+s1_range(1)+s2_range(1)];
ax = -cropSize(1)/2+1:cropSize(1)/2;

figure, imagesc(ax,ax,recos_crop(:,:,1,2),gsrange)
axis equal tight, colormap gray
hold on
rectangle('Position',[cent(1)-rad(1),cent(2)-rad(1),rad(1)*2,rad(1)*2],'Curvature',[1,1],...
    'EdgeColor','r','LineWidth',2)
rectangle('Position',[cent(1)-rad(2),cent(2)-rad(2),rad(2)*2,rad(2)*2],'Curvature',[1,1],...
    'EdgeColor','r','LineWidth',2)
rectangle('Position',[cent(1)-rad(3),cent(2)-rad(3),rad(3)*2,rad(3)*2],'Curvature',[1,1],...
    'EdgeColor','r','LineWidth',2)

%% Processing projections pass 1
% first pass:
%   - flat and dark correction
%   - stitching in x
%   - generating mean projections
%       - can generate a ring correction, however in the present
%       implementation of ProjectionProcessing_pass2, this is not used and
%       a different ring correction is made
% Note: current implementation assumes all hs have same x- stitch positions
%   it would be fairly easy input a matrix of values as well
manstitchposx = [
        207.6    1848.3    1844.6    1848.3;  % 1
        205.9    1848.3    1848.3    1848.3;  % 2
        205.3    1848.6    1848.6    1847.7;  % 3
        205.4    1848.1    1848.6    1847.8;  % 4
        205.0    1848.3    1848.3    1848.3;  % 5
        205.0    1848.5    1846.0    1847.6;  % 6
        203.6    1848.3    1848.3    1848.3;  % 7
        204.6    1848.3    1848.3    1848.3;  % 8
    ];
projsavedir = ProjectionProcessing_pass1(paramfile,manstitchposx);

%% Automatically find height step stitching positions (all height steps)
% note: this is much more robust when using full stitched projections
%
% currently, OverlapFinderY runs after ProjectionProcessing_pass1
% you can easily modify it to run before, look at lines imr1 = imread()
% and replace with h5read of projection (look to OverlapFinderX for snippet
% to use)
% 
% Alternative would be to re-write to process the stitched projections as
% you go in the y overlap finder
[writedir] = OverlapFinderY(paramfile);

%% Check height step stitching positions in some projections

manstitchposy = [1835,1843,1843,1847,1847,1848,1839];
proj_nr = 1:500:4001;

projsavedir = y_stitching_check_proj(paramfile,manstitchposy,proj_nr);

%% Check height step stitching positions manually
% % Requires ProjectionProcessing_pass1 has been run for hs in question
% Note: this section does not quite work yet, it should be adjusted

% % set up a directory for tests
testdir = [projdir samplename filesep 'stitchpos_tests' filesep];
if not(isfolder(testdir)); mkdir(testdir); end

% % Which overlap?
hs = 7; % 1 is 1-2, 2 is 2-3, etc.§

% stitch position for this overlap
manstitchposy = 1839;
ycrop1 = 1941:1956; % region in overlap, i.e. larger than manstitchposy (length 16)
xcrop = 6000; % number of pixels to crop projections from both sides
% projsavedir = ProjectionProcessing_pass2_YCheck(paramfile,manstitchposy,hs,ycrop1,ycrop2);
projsavedir = ProjectionProcessing_pass2_YCheck_v2(paramfile,manstitchposy,hs,ycrop1,xcrop);

tmp1 = dir([projsavedir 'proj1_*.tif']);
tmp2 = dir([projsavedir 'proj2_*.tif']);
tmpim = imread([tmp1(1).folder filesep tmp1(1).name]);
[sy,sx] = size(tmpim);
vol1 = zeros(sy,sx,length(tmp1),'single');
vol2 = zeros(sy,sx,length(tmp2),'single');
for i = 1:length(tmp1)
    vol1(:,:,i) = imread([tmp1(i).folder filesep tmp1(i).name]);
end
for i = 1:length(tmp2)
    vol2(:,:,i) = imread([tmp2(i).folder filesep tmp2(i).name]);
end

reco1 = zeros(sy, sx, sx, 'single');
reco2 = zeros(sy, sx, sx, 'single');

parfor y = 1:sy
    reco1(y,:,:) = SingleGridrecReconstruction(...
        [testdir 'overlap1_slice' num2str(y)],...
        squeeze(vol1(y,:,:)),...
        angles,...
        pixsize_mm,...
        pythonscript_fullpath);
    reco2(y,:,:) = SingleGridrecReconstruction(...
        [testdir 'overlap2_slice' num2str(y)],...
        squeeze(vol2(y,:,:)),...
        angles,...
        pixsize_mm,...
        pythonscript_fullpath);
end

fprintf(['writing volume 1 to ' testdir 'overlap1_reco.mhd\n'])
mha_write(permute(reco1, [3 2 1]), [testdir 'overlap1_reco'])

fprintf(['writing volume 2 to ' testdir 'overlap2_reco.mhd\n'])
mha_write(permute(reco2, [3 2 1]), [testdir 'overlap2_reco'])

fprintf('deleting temporaries\n')
for h = 1:2
    for y = 1:sy
        delete([testdir 'overlap' num2str(h) '_slice' num2str(y) '.h5']);
        delete([testdir 'overlap' num2str(h) '_slice' num2str(y) ...
            '_reco.h5']);
    end
end

%% Check translation of rotation center
% We observed a displacement in the two subsequent height steps.
% The displacement is determined with elastix.
% Here the projections are individually shifted, to result in a translation
% of the rotation center.
% Previous step needs to be run first and y-stitching should already
% have been adjusted.

% tx and ty as output by elastix for EulerTransform,
% i.e. for the mapping fixed->moving image
% (TransformParameters thetax, thetay, thetaz, tx, ty, tz)
tx = 8.25;
ty = -5.29;

angles_rad = angles * pi / 180.
d_alpha = -ty * cos(angles_rad) + tx * cos(angles_rad);

vol2_shift = zeros(sy, sx, length(angles), 'single');
parfor i = 1:length(angles)
    vol2_shift(:,:,i) = subpixelshift(vol2(:,:,i), 0, d_alpha(i));
end

reco2_shift = zeros(sy, sx, sx, 'single');
parfor y = 1:sy
    reco2_shift(y,:,:) = SingleGridrecReconstruction(...
        [testdir 'overlap2_shift_slice' num2str(y)],...
        squeeze(vol2_shift(y,:,:)),...
        angles, pixsize_mm,...
        pythonscript_fullpath);
end

mha_write(permute(reco2_shift, [3,2,1]), [testdir 'overlap2_shift_reco']);

%% Test filtering and ring correction parameters
% % set up a directory for tests
testdir = [projdir samplename filesep 'db_ringcorr_tests' filesep];
if not(isfolder(testdir)); mkdir(testdir); end

% % height and crop to use
this_hs = 3;
this_ycrop = 1024-7:1024+8;

cor = 205;
s1x = 1848.6;
s2x = 1848.6;
s3x = 1847.7;
manstitchposx = [cor,s1x,s2x,s3x];
projdir = ProjectionProcessingManualEntry(paramfile,this_hs,manstitchposx,...
        this_ycrop,0,'none');
[sinoblock] = LoadCroppedSinoblock(projdir);

% % testing paganin filtering
test_photonEnergy = photonEnergy; % [keV]
test_det_dist_mm = det_dist_mm; % [mm]
test_pixsize_mm = pixsize_mm; % [mm]
test_lambda_mm = lambda_from_E(test_photonEnergy*1e3)*1e3; % [mm] (note: lambda_from_E takes [eV] and gives [m])

cropSize = [14000,14000];

dbrange = [0,50,100,200,500,1000];
recos_crop = zeros(cropSize(2),cropSize(1),length(dbrange),'single');
for db = 1:length(dbrange)
    thisdb = dbrange(db);
    if thisdb == 0
    thissinoblock = sinoblock;
    else
    thissinoblock = paganin_filter_stack(sinoblock,...
        thisdb,test_pixsize_mm,test_lambda_mm,test_det_dist_mm);
    end
    fullsino = squeeze(thissinoblock(round(end/2),:,:));
    fullsino = fullsino(5:end-5,:);
    reco = SingleGridrecReconstruction([testdir ...
        'reco_db_' num2str(thisdb)],fullsino,angles,test_pixsize_mm,...
        pythonscript_fullpath);
    recos_crop(:,:,db) = cropaboutcenter(reco,cropSize);
end

gsrange = [-0.01,0.05];
figure, imshow3D(squeeze(recos_crop),gsrange)

% % testing ring correction
% baseline, no correction
thissinoblock = paganin_filter_stack(sinoblock,...
    filterwidth,pixsize_mm,lambda_mm,det_dist_mm);
fullsino = squeeze(thissinoblock(round(end/2),:,:));
fullsino = fullsino(5:end-5,:);
reco_nrc = SingleGridrecReconstruction([testdir ...
    'reco_nrc_db_' num2str(filterwidth)],fullsino,angles,pixsize_mm,...
    pythonscript_fullpath);


% this is the currently used ring correction:
mproj = mean(sinoblock,3);
rproj = single(mproj-imgaussfilt(mproj,50,'Padding','symmetric'));
rproj(abs(rproj)>0.1) = 0;

thissinoblock = paganin_filter_stack(sinoblock-rproj,...
    filterwidth,pixsize_mm,lambda_mm,det_dist_mm);
fullsino = squeeze(thissinoblock(round(end/2),:,:));
fullsino = fullsino(5:end-5,:);
reco_rc = SingleGridrecReconstruction([testdir ...
    'reco_rc_db_' num2str(filterwidth)],fullsino,angles,pixsize_mm,...
    pythonscript_fullpath);

h5create([projsavedir 'angles.h5'], '/angles', size(angles), 'Datatype', 'single');
h5write([projsavedir 'angles.h5'], '/angles', single(angles));

%% Projection processing
% second pass:
%   - stitching in y
%   - (optional) ring correction -- in this case lines 92-97
%   - (optional) filtering

manstitchposy = [1835,1843,1843,1847,1847,1848,1839];
translation = [-0.664023, -5.645471,  1.051685;...
    0.236059, 0.497671,   -3.858646;...
    0.528731, 0.774861,   -3.760332;...
    -0.255347, 1.581965,   -2.778743;...
    3.836355, 4.420312,   -5.637497;...
    7.248668, 7.995099,   -2.645331;...
    -5.286473, 8.250312,   -1.135933];
projsavedir = ProjectionProcessing_pass2(paramfile,manstitchposy,...
    translation);

%% Decide output scaling and cropping
% would recommend loading sinograms from various positions, reconstructing
% them and using that to determine what to use for "outputcrop" and "outputgrayscale" 

%% Reconstruction


%% Binning




