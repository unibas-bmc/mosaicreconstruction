function [proj,figdir] = PreviewStitchProjectionXY(paramfile,projNo,medkernel,gausskernel)
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

h5ImagePath = infoStruct.h5ImagePath;
h5AnglePath = infoStruct.h5AnglePath;
%% 0.4 Set up directories
procdir = [projdir samplename filesep];
figdir = [procdir filesep 'figs' filesep];
if not(isfolder(figdir)); mkdir(figdir); end
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
%% 1.0 Make stitched projection
fprintf('Loading individual projections...\n')
nflats = 4; fstep = 1;

proj = zeros(datsize(1),datsize(2),nhs,nrings);
proj_mf = zeros(datsize(1),datsize(2),nhs,nrings);
proj_mfgf = zeros(datsize(1),datsize(2),nhs,nrings);
projp180 = zeros(datsize(1),datsize(2),nhs,nrings);
projp180_mf = zeros(datsize(1),datsize(2),nhs,nrings);
projp180_mfgf = zeros(datsize(1),datsize(2),nhs,nrings);

for h = 1:nhs
    for r = 1:nrings 
        %% 1.1 load flats/darks
        tmp = rot90(h5read([rawbasedir ringnames{h,r} filesep 'pre_ref.nxs'],'/ref_0',...
            [1,1,1],[datsize(1),datsize(2),nflats],[1,1,fstep]));
        tmp1 = rot90(h5read([rawbasedir ringnames{h,r} filesep 'post_ref.nxs'],'/ref_0',...
            [1,1,1],[datsize(1),datsize(2),nflats],[1,1,fstep]));
        flats = double(cat(3,tmp,tmp1)); clear tmp tmp1
        mflats = median(flats,3);
        darks = rot90(h5read([rawbasedir ringnames{h,r} filesep 'post_dark.nxs'],'/dark_0'));
        mddarks = double(median(darks,3));
        
        %% 1.2 Loading projections
        tmp = double(rot90(h5read([rawbasedir ringnames{h,r} filesep ringnames{h,r} '.nxs'],h5ImagePath,...
            [1,1,projNo],[datsize(1),datsize(2),1])));
        tmp = flipud((tmp-mddarks)./(mflats-mddarks));
        proj(:,:,h,r) = tmp;
        proj_mf(:,:,h,r) = medfilt2(tmp,medkernel);
        proj_mfgf(:,:,h,r) = imgaussfilt(medfilt2(tmp,medkernel),gausskernel);
        tmp = double(rot90(h5read([rawbasedir ringnames{h,r} filesep ringnames{h,r} '.nxs'],h5ImagePath,...
                [1,1,ip180+projNo-1],[datsize(1),datsize(2),1])));
        tmp = flipud((tmp-mddarks)./(mflats-mddarks));
        projp180(:,:,h,r) = tmp;
        projp180_mf(:,:,h,r) = medfilt2(tmp,medkernel);
        projp180_mfgf(:,:,h,r) = imgaussfilt(medfilt2(tmp,medkernel),gausskernel);

        
        
    end
end
%% 1.3 Guess stitch positions from motor postitions
ystitch_guess = abs((tzpos_all(2:end,:)-tzpos_all(1:end-1,:))/pixsize_mm);
xstitch_guess = zeros(nhs,nrings-1);
for iii = 1:nhs
    for ii = 1:nrings-1
        this_xshift = (txpos_all(iii,ii)-txpos_all(iii,ii+1));
        xstitch_guess(iii,ii) = this_xshift/pixsize_mm;   
    end
end
cor_guess = zeros(1,nhs);
for iii = 1:nhs
    this_cor = txpos_all(iii,1);
    cor_guess(iii) = 2*round((datsize(2)/2)+(this_cor/pixsize_mm));
end
ystitch = 0;
if ~isempty(ystitch_guess)
    ystitch = round(mean(ystitch_guess,'all'));
end
xstitch = round(mean(xstitch_guess,'all'));
cor = round(mean(cor_guess));
%% 1.4 Build x-projections
% build x projections
stitchposx = [0,cumsum(xstitch*ones(1,nrings-1))];
halfwidth = ceil(datsize(2)+sum(xstitch*ones(1,nrings-1))); % width before stitching off-axis
fullwidth = halfwidth*2-ceil(cor);
osy = datsize(1);
if 1 == 1
halfmask = zeros(osy,halfwidth,nrings);
for r = 1:nrings
    if r == 1
        halfmask(:,floor(stitchposx(r))+1:floor(stitchposx(r))+datsize(2)-1,r) = ones(osy,datsize(2)-1);
        tmp = floor(stitchposx(r+1))+1:floor(stitchposx(r))+datsize(2)-1;
        halfmask(:,tmp,r) = ones(osy,length(tmp)).*linspace(1,0,length(tmp));
    elseif r == nrings
        halfmask(:,floor(stitchposx(r))+2:floor(stitchposx(r))+datsize(2),r) = ones(osy,datsize(2)-1);
        tmp = floor(stitchposx(r))+1:floor(stitchposx(r-1))+datsize(2)-1;
        halfmask(:,tmp,r) = ones(osy,length(tmp)).*linspace(0,1,length(tmp));
    else
        halfmask(:,floor(stitchposx(r))+2:floor(stitchposx(r))+datsize(2)-1,r) = ones(osy,datsize(2)-2);
        tmp = floor(stitchposx(r+1))+1:floor(stitchposx(r))+datsize(2)-1;
        halfmask(:,tmp,r) = ones(osy,length(tmp)).*linspace(1,0,length(tmp));
        tmp = floor(stitchposx(r))+1:floor(stitchposx(r-1))+datsize(2)-1;
        halfmask(:,tmp,r) = ones(osy,length(tmp)).*linspace(0,1,length(tmp));
    
    end
end
% build a mask for the 0 plus 180 degrees projection
fullmask = zeros(osy,fullwidth,2);
fullmask(:,1:halfwidth-1,1) = ones(osy,halfwidth-1);
fullmask(:,halfwidth-ceil(cor)+2:end,2) = ones(osy,halfwidth-1);
tmp = halfwidth-ceil(cor)+2:halfwidth;
fullmask(:,halfwidth-ceil(cor)+2:halfwidth,1) = ones(osy,length(tmp)).*linspace(1,0,length(tmp));
fullmask(:,halfwidth-ceil(cor)+2:halfwidth,2) = ones(osy,length(tmp)).*linspace(0,1,length(tmp));
end
ringproj = zeros(datsize(1),fullwidth,nhs);
for h = 1:nhs
    rim0 = squeeze(proj_mfgf(:,:,h,:));
    rim180 = squeeze(projp180_mfgf(:,:,h,:));
    
    ebump = 10;
    for i = fliplr(1:nrings-1)
        olregi = floor(stitchposx(i+1))+1:floor(stitchposx(i))+datsize(2);
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
        proj(:,floor(stitchposx(i))+1:floor(stitchposx(i))+datsize(2),i) = rim0(:,:,i);
        projp180(:,floor(stitchposx(i))+1:floor(stitchposx(i))+datsize(2),i) = rim180(:,:,i);
    end
    proj = BlendProjection(proj,halfmask);
    halfprojp180 = BlendProjection(projp180,halfmask);
    proj = fliplr(proj);
    morl = mean(proj(:,halfwidth-ceil(cor)+1:end),'all');
    morr = mean(halfprojp180(:,1:ceil(cor)),'all');
    proj = (morr/morl)*proj;
    
    tmpproj = zeros(osy,fullwidth,2);
    tmpproj(:,1:halfwidth,1) = proj;
    tmpproj(:,halfwidth-ceil(cor)+1:end,2) = halfprojp180;
    fullproj = (BlendProjection(tmpproj,fullmask));
    
    % filter stitched projection
    bump = 5;
    fullproj = fullproj(:,bump+1:end-bump); % knock off edges
    fullproj = padarray(fullproj,[0,bump],'symmetric','both');
    ringproj(:,:,h) = flipud(fullproj);
end
%% 1.5 Build y-projection
% Note: "halfheight" is used because this is adapted from the
% x-stitching ("halfwidth"), but there is no "other half" in y-stitching, 
% sorry it isn't the best naming convention
stitchposy = [0,cumsum(ystitch*ones(1,nhs-1))];
halfheight = ceil(datsize(2)+sum(ystitch*ones(1,nhs-1))); % width before stitching off-axis
if 1 == 1
halfmask = zeros(halfheight,fullwidth,nhs);
for r = 1:nhs
    if r == 1
        halfmask(floor(stitchposy(r))+1:floor(stitchposy(r))+datsize(2)-1,:,r) = ones(datsize(2)-1,fullwidth);
        tmp = floor(stitchposy(r+1))+1:floor(stitchposy(r))+datsize(2)-1;
        halfmask(tmp,:,r) = ones(length(tmp),fullwidth).*linspace(1,0,length(tmp))';
    elseif r == nhs
        halfmask(floor(stitchposy(r))+2:floor(stitchposy(r))+datsize(2),:,r) = ones(datsize(2)-1,fullwidth);
        tmp = floor(stitchposy(r))+1:floor(stitchposy(r-1))+datsize(2)-1;
        halfmask(tmp,:,r) = ones(length(tmp),fullwidth).*linspace(0,1,length(tmp))';
    else
        halfmask(floor(stitchposy(r))+2:floor(stitchposy(r))+datsize(2)-1,:,r) = ones(datsize(2)-2,fullwidth);
        tmp = floor(stitchposy(r+1))+1:floor(stitchposy(r))+datsize(2)-1;
        halfmask(tmp,:,r) = ones(length(tmp),fullwidth).*linspace(1,0,length(tmp))';
        tmp = floor(stitchposy(r))+1:floor(stitchposy(r-1))+datsize(2)-1;
        halfmask(tmp,:,r) = ones(length(tmp),fullwidth).*linspace(0,1,length(tmp))';
    
    end
end
end

rim0 = ringproj;
ebump = 10;
for i = fliplr(1:nhs-1)
    olregi = floor(stitchposy(i+1))+1:floor(stitchposy(i))+datsize(2);
    olregi = datsize(2)-length(olregi)+1:datsize(2);
    olregi = olregi(ebump):olregi(end-ebump);
    olregip1 = (1:length(olregi))+ebump;

    mor0i = mean(rim0(olregi,:,i),'all');
    mor0ip1 = mean(rim0(olregip1,:,i+1),'all');
        
    rim0(:,:,i) = (mor0ip1/mor0i)*rim0(:,:,i);
end
proj = zeros(halfheight,fullwidth,nhs);
for i = 1:nhs
    proj(floor(stitchposy(i))+1:floor(stitchposy(i))+datsize(2),:,i) = rim0(:,:,i);
end
proj = BlendProjection(proj,halfmask);
proj = flipud(proj);
%% 2.0 Save some figures
prepimage = @(im,vr) uint8(255*(double(im)-vr(1))/(vr(2)-vr(1)));

proj_b10 = imresize(proj,0.1);

sty = size(proj,1)-stitchposy(2:end);
stx = [size(proj,2)/2-fliplr(stitchposx),size(proj,2)/2+stitchposx];
stx = unique(stx);

figure, imagesc(1:size(proj,2),1:size(proj,1),proj,[0.75,1.0]), axis equal tight, colormap(flipud(bone))
hold on
plot(repmat(stx,[2,1]),repmat([1,size(proj,1)],[length(stx),1])','r-')
plot(repmat([1,size(proj,2)],[length(sty),1])',repmat(sty,[2,1]),'r-')
f = gcf; f.Position = [137.6667  161.6667  424.0000  385.3333];
saveas(f,[figdir 'proj' num2str(projNo) '_axison_stitched_stitchpos_75to100.pdf'])

figure, imagesc(1:size(proj,2),1:size(proj,1),proj,[0.75,0.85]), axis equal tight, colormap(flipud(bone))
hold on
plot(repmat(stx,[2,1]),repmat([1,size(proj,1)],[length(stx),1])','r-')
plot(repmat([1,size(proj,2)],[length(sty),1])',repmat(sty,[2,1]),'r-')
f = gcf; f.Position = [137.6667  161.6667  424.0000  385.3333];
saveas(f,[figdir 'proj' num2str(projNo) '_axison_stitched_stitchpos_75to85.pdf'])

figure, imagesc(1:size(proj,2),1:size(proj,1),proj_b10,[0.75,1.0]), axis equal tight, colormap(flipud(bone))
hold on
plot(repmat(stx,[2,1]),repmat([1,size(proj,1)],[length(stx),1])','r-')
plot(repmat([1,size(proj,2)],[length(sty),1])',repmat(sty,[2,1]),'r-')
axis off
f = gcf; f.Position = [137.6667  161.6667  424.0000  385.3333];
saveas(f,[figdir 'proj' num2str(projNo) '_stitched_stitchpos_75to100.pdf'])

figure, imagesc(1:size(proj,2),1:size(proj,1),proj_b10,[0.75,0.85]), axis equal tight, colormap(flipud(bone))
hold on
plot(repmat(stx,[2,1]),repmat([1,size(proj,1)],[length(stx),1])','r-')
plot(repmat([1,size(proj,2)],[length(sty),1])',repmat(sty,[2,1]),'r-')
axis off
f = gcf; f.Position = [137.6667  161.6667  424.0000  385.3333];
saveas(f,[figdir 'proj' num2str(projNo) '_stitched_stitchpos_75to85.pdf'])

imwrite(ind2rgb(imcomplement(prepimage(proj_b10,[0.75,1.0])),bone),[figdir 'proj' num2str(projNo) '_stitched_75to100_invbone.tif']);
imwrite(ind2rgb(imcomplement(prepimage(proj_b10,[0.75,0.85])),bone),[figdir 'proj' num2str(projNo) '_stitched_75to85_invbone.tif']);
imwrite(imcomplement(prepimage(proj_b10,[0.75,1.0])),[figdir 'proj' num2str(projNo) '_stitched_75to100_inv.tif']);
imwrite(imcomplement(prepimage(proj_b10,[0.75,0.85])),[figdir 'proj' num2str(projNo) '_stitched_75to85_inv.tif']);
imwrite(prepimage(proj_b10,[0.75,1.0]),[figdir 'proj' num2str(projNo) '_stitched_75to100.tif']);
imwrite(prepimage(proj_b10,[0.75,0.85]),[figdir 'proj' num2str(projNo) '_stitched_75to85.tif']);
end

