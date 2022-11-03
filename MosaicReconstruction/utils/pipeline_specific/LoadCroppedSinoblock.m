function [sinoblock] = LoadCroppedSinoblock(projdir)
dirinfo = dir([projdir 'proj*.tif']);

sz = length(dirinfo);
t = Tiff([projdir dirinfo(1).name]);
sx = t.getTag('ImageWidth');
sy = t.getTag('ImageLength');
nblocks = numberOfStrips(t);
blocksize = t.getTag('RowsPerStrip');

sinoblock = zeros(sy,sx,sz,'single');
for i = 1:sz
    t = Tiff([projdir dirinfo(i).name]);
    %sinoblock(:,:,i) = readEncodedStrip(t,1);
    sinoblock(:,:,i) = read(t);
end

end

