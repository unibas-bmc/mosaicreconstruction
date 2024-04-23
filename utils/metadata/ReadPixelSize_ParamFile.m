function [pixsize_mm] = ReadPixelSize_ParamFile(paramfile)

fid = fopen(paramfile);
infoStruct = textscan(fid, '%s %s','Delimiter','\t','CommentStyle','//');
infoStruct = cell2struct(infoStruct{2},infoStruct{1},size(infoStruct,1));
samplename = infoStruct.samplename;
tmp = split(infoStruct.heightstep,'-');
hs_range = str2double(tmp{1}):str2double(tmp{2});
nhs = length(hs_range);
infofile = infoStruct.infofile;
rawdatapath = infoStruct.rawdatapath;
fclose(fid);

T = readtable(infofile); 
scandir = T.scandir{1};
this_suffix = T.suffix{1};
if this_suffix == 0
    scansuffix = '';
else
    scansuffix = ['_' num2str(this_suffix)];
end

parfile = [rawdatapath scandir scansuffix filesep scandir scansuffix '.par'];


myStr1 = 'IMAGE_PIXEL_SIZE_1 = '; lenStr1 = length(myStr1);

fid=fopen(parfile);
i = 1;
tline = fgetl(fid);
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    if contains(tline,myStr1)
        break
    end
end
fclose(fid);

tline = tline(lenStr1+1:end);
newStr = split(tline,' ');
pixsize_um = str2double(newStr(1));
pixsize_mm = pixsize_um*1e-3;
end

