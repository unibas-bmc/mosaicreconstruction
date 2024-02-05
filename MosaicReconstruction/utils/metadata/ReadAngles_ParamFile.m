function [angles,ip180] = ReadAngles_ParamFile(paramfile,h5AnglePath)


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

nxsfile = [rawdatapath scandir scansuffix filesep scandir scansuffix '.nxs'];


find_index = @(v,val) find(abs(v-val) == min(abs(v-val)),1);

angles = h5read(nxsfile,h5AnglePath);
angles = angles - angles(1);
if angles(2) > 0
    ip360 = find_index(angles,360);% index of projection with angle 360 degrees
else
    ip360 = find_index(angles,-360);% index of projection with angle -360 degrees
end
ip180 = floor(ip360/2); % index of projection with angle 180 degrees
angles = angles(1:ip180);

end
