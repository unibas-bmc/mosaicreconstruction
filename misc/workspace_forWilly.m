%% Paths and functions
addpath(genpath('../reco_scripts'))

cropaboutcenter = @(im,roisize) im(floor(size(im,1)/2)-floor(roisize(2)/2)+1:floor(size(im,1)/2)+floor(roisize(2)/2),...
    floor(size(im,2)/2-roisize(1)/2)+1:floor(size(im,2)/2+roisize(1)/2));
cropto = @(im,roi) im(roi(3):roi(4),roi(1):roi(2));
cropto3d = @(im,roi) im(roi(3):roi(4),roi(1):roi(2),roi(5):roi(6));
croptostack = @(im,roi) im(roi(3):roi(4),roi(1):roi(2),:);
prepimage = @(im,vr) uint8(255*(double(im)-vr(1))/(vr(2)-vr(1)));
prepimagei16 = @(im,vr) uint16((2^16-1)*(double(im)-vr(1))/(vr(2)-vr(1)));

%% Manual overlap finder
paramfile = ['/media/griffin/_home1/mosaic_reconstruction/mousebrain/' ...
    'param_files/mouse4_perf_eth_hs5_param.txt'];

procdir = overlap_finder(paramfile);

readdir = projection_processing_manualoverlap(paramfile);

[projvol,mprojvol] = read_projections_manualoverlap(paramfile);

% % apply any desired paganin filtering and ring correction
photon_energy = 27e3; % [eV]
det_dist_mm = 50;
lambda = lambda_from_E(photon_energy); % [m]
lambda_mm = lambda*1e3;
deltabeta = 48; % delta/beta

projvol_pag = projvol;
for i1 = 1:size(projvol,4)
    projvol_pag(:,:,:,i1) = paganin_filter_stack(projvol(:,:,:,i1),...
        deltabeta,pixsize_mm,lambda_mm,det_dist_mm);
end

% % Load overlap finder results
tmpdir = [procdir 'parameters' filesep];
load([tmpdir 'stitching_parameters.mat'],'ol_final_subpix','corpix_subpix','fullwidth','stitchpos','halfwidth')

% % read motor positions
motorPositions = read_motors(paramfile);
angles = motorPositions.angles;
nrings = motorPositions.nrings;
ip180 = length(angles);
pixsize_mm = motorPositions.pixsize_mm;

% % set up a directory for tests
testdir = [procdir 'test_recos' filesep];
if not(isfolder(testdir)); mkdir(testdir); end

% % check center of rotation
sliceNo = 9;
corRange = 200:0.5:204;
% note: motor position would be cor_guess, found position would be
% cor_subpix

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
    
    fullsino = (blendproj(fullsino,fullmask))';
    fullsino = padarray(fullsino,[padSize,0],'symmetric','both'); % padding since this is local

    
    reco = single_gridrec_reconstruction([testdir ...
        'localreco_COR_' num2str(this_cor)],fullsino,angles,pixsize_mm);
    reco = reco(padSize+1:end-padSize,padSize+1:end-padSize);

    recos_cor(:,:,i1) = cropaboutcenter(reco,cropSize);
end

figure, imshow3D(recos_cor,prctile(recos_cor,[3,97],'all')')

% % check stitching positions of each ring
% % Ring 1
this_cor = 202.6;

olpix1_range = 1842.6:1852.6;
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

    halfproj = blendproj(proj,halfmask); % note 10.03.2021: blendproj should be improved
    halfprojp180 = blendproj(projp180,halfmask);
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
    fullproj = (blendproj(tmpproj,fullmask))';
    fullproj = fullproj(2:end-1,:);
    
    fullproj = padarray(fullproj,[padSize,0],'symmetric','both'); % padding since this is local

    reco = single_gridrec_reconstruction([testdir ...
        'localreco_ol_r1_' num2str(olpix1_range(i1))],fullproj,angles,pixsize_mm);
    reco = reco(padSize+1:end-padSize,padSize+1:end-padSize);

    recos_crop(:,:,i1) = cropaboutcenter(reco,cropSize);
end
figure, imshow3D(recos_crop,prctile(recos_crop,[3,97],'all')')

cent = [0,0];
rad = 2048-this_cor/2;
ax = -cropSize(1)/2+1:cropSize(1)/2;
figure, imagesc(ax,ax,recos_crop(:,:,round(end/2)),prctile(recos_crop,[3,97],'all')')
axis equal tight, colormap gray
hold on
%viscircles(cent,rad,'EnhanceVisibility',false);
rectangle('Position',[cent(1)-rad,cent(2)-rad,rad*2,rad*2],'Curvature',[1,1],...
    'EdgeColor','r')

c1 = [3750+1000,3750+2000,3750-500,3750+500];
c1p5 = [3750+1500,3750+2500,3750-500,3750+500];
c2 = [3750+2000,3750+3000,3750-500,3750+500];
recos_c1 = croptostack(recos_crop,c1);
recos_c1p5 = croptostack(recos_crop,c1p5);
recos_c2 = croptostack(recos_crop,c2);

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.01 0.05], [0.05 0.025]);

for i1 = 1:length(olpix1_range)
figure, subplot(1,2,1)
imagesc(ax(c1(1):c1(2)),ax(c1(3):c1(4)),...
    recos_c1(:,:,i1),prctile(recos_c1(:,:,i1),[1,99],'all')')
axis equal tight, colormap gray
title('ring 1')
subplot(1,2,2)
imagesc(ax(c2(1):c2(2)),ax(c2(3):c2(4)),...
    recos_c2(:,:,i1),prctile(recos_c2(:,:,i1),[1,99],'all')')
axis equal tight, colormap gray
title('ring 2')
sgtitle(['Stitch Position = ' num2str(olpix1_range(i1)) ' pixels'])
end

figure, imshow3D(recos_c2,prctile(recos_c2,[1,99],'all')')
title('ring 2')
figure, imshow3D(recos_c1,prctile(recos_c1,[1,99],'all')')
title('ring 1')
figure, imshow3D(recos_c1p5,prctile(recos_c1p5,[1,99],'all')')
title('border')


% % Tweak any
cor_range = 201.6:204.6;
s1_range = 1843.6:2:1851.6;
s2_range = 1847.9;
s3_range = 1847.9;

min_size = ceil(2048+max([0,cumsum([min(s1_range),min(s2_range),min(s2_range)])]))*2-ceil(min(cor_range));

padSize = 0;
cropSize = [14000,14000];
recos_crop = zeros(cropSize(2),cropSize(1),length(cor_range),length(s1_range),...
    length(s2_range),length(s3_range),'single');
for i1 = 1:length(cor_range)
for i2 = 1:length(s1_range)
for i3 = 1:length(s2_range)
for i4 = 1:length(s3_range)
    these_olpix = [s1_range(i2),s2_range(i3),s2_range(i4)];
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

    halfproj = blendproj(proj,halfmask); % note 10.03.2021: blendproj should be improved
    halfprojp180 = blendproj(projp180,halfmask);
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
    fullproj = (blendproj(tmpproj,fullmask))';
    fullproj = fullproj(2:end-1,:);
    
    fullproj = padarray(fullproj,[padSize,0],'symmetric','both'); % padding since this is local

    reco = single_gridrec_reconstruction([testdir ...
        'cor_' num2str(cor_range(i1)) '_s1_' num2str(s1_range(i2)) ...
        '_s2_' num2str(s2_range(i3)) '_s3_' num2str(s3_range(i4)) ],...
        fullproj,angles,pixsize_mm);
    reco = reco(padSize+1:end-padSize,padSize+1:end-padSize);

    recos_crop(:,:,i1,i2,i3,i4) = cropaboutcenter(reco,cropSize);
end
end
end
end

cent = [0,0];
rad = [2048-cor_range(1)/2,2048-cor_range(1)/2+s1_range(1),...
    2048-cor_range(1)/2+s1_range(2)+s2_range(1)];
ax = -cropSize(1)/2+1:cropSize(1)/2;

figure, imagesc(ax,ax,recos_crop(:,:,1,2),prctile(recos_crop,[3,97],'all')')
axis equal tight, colormap gray
hold on
rectangle('Position',[cent(1)-rad(1),cent(2)-rad(1),rad(1)*2,rad(1)*2],'Curvature',[1,1],...
    'EdgeColor','r')
rectangle('Position',[cent(1)-rad(2),cent(2)-rad(2),rad(2)*2,rad(2)*2],'Curvature',[1,1],...
    'EdgeColor','r')
rectangle('Position',[cent(1)-rad(3),cent(2)-rad(3),rad(3)*2,rad(3)*2],'Curvature',[1,1],...
    'EdgeColor','r')


c1 = [7000+850,7000+1350,7000+850,7000+1350];
c2 = [7000+2000,7000+2500,7000+2000,7000+2500];
c3 = [7000+3200,7000+3700,7000+3200,7000+3700];
c4 = [7000+4050,7000+4550,7000+4050,7000+4550];
tmpfunc = @(im,roi) im(roi(3):roi(4),roi(1):roi(2),:,:);
recos_c1 = tmpfunc(recos_crop,c1);
recos_c2 = tmpfunc(recos_crop,c2);
recos_c3 = tmpfunc(recos_crop,c3);
recos_c4 = tmpfunc(recos_crop,c4);

figure, imshow3D(squeeze(recos_c1(:,:,2,:)),prctile(recos_c1,[1,99],'all')')
title('ring 1')
figure, imshow3D(squeeze(recos_c2(:,:,2,:)),prctile(recos_c2,[1,99],'all')')
title('ring 2')
figure, imshow3D(squeeze(recos_c3(:,:,2,:)),prctile(recos_c3,[1,99],'all')')
title('ring 3')
figure, imshow3D(squeeze(recos_c4(:,:,2,:)),prctile(recos_c4,[1,99],'all')')
title('ring 4')

c1_1 = [61,160,101,200];
c2_1 = [216,315,131,230];
c3_1 = [116,215,196,296];
c4_1 = [371,470,71,170];
recos_c1_1 = tmpfunc(recos_c1,c1_1);
recos_c2_1 = tmpfunc(recos_c2,c2_1);
recos_c3_1 = tmpfunc(recos_c3,c3_1);
recos_c4_1 = tmpfunc(recos_c4,c4_1);

figure, imshow3D(squeeze(recos_c1_1(:,:,2,:)),prctile(recos_c1,[1,99],'all')')
title('ring 1')
figure, imshow3D(squeeze(recos_c2_1(:,:,2,:)),prctile(recos_c2,[1,99],'all')')
title('ring 2')
figure, imshow3D(squeeze(recos_c3_1(:,:,2,:)),prctile(recos_c3,[1,99],'all')')
title('ring 3')
figure, imshow3D(squeeze(recos_c4_1(:,:,2,:)),prctile(recos_c4,[1,99],'all')')
title('ring 4')

