function [datsize] = ReadScanSize(filename,h5path)
tmp = h5info(filename,h5path);
datsize = tmp.Dataspace.Size; % size of images (usually 2048x2048)
end

