function[XYZ,CHECK] = XYZexporter(GRIDDATA,INFO)
%take matrix format data and exports to XYZ format.
%imports are GRIDDATA: the 2-D matrix of data to be exported and
%INFO which is a struct field with elements:
%xmin, xmax, ymin, ymax, xstep, ystep, xlabel, ylabel.
%Each field is a scalar with the exception of xlabel and ylabel, which are
%vectors. CHECK is a return value that gives the state of the function
%   Detailed explanation goes here

%initialize variables
CHECK=0;
idx=1;
L=sum(sum(~isnan(GRIDDATA),2));
XYZ=zeros(L,3);
INFO.ylabel=fliplr(INFO.ylabel);
CHECK=1;

%iterate Y
for i=1:length(GRIDDATA(:,1))
    %iterate X
    for j=1:length(GRIDDATA(1,:))
       if isnan(GRIDDATA(i,j)) == 0
     XYZ(idx,1)=INFO.xlabel(j);
     XYZ(idx,2)=INFO.ylabel(i);
     XYZ(idx,3)=GRIDDATA(i,j);
     idx=idx+1;
       end
    end
end
CHECK=2;
fname=sprintf('Grid_Export_%s.dat',datetime('now','format','yyyy-MM-dd_HH-mm-ss'));
[file,path]=uiputfile(fname,'Select File Name');
%add .dat to the end of the file, if it's not already there
if strcmp(file(end-3:end),'.dat') == 0
    file=strcat(file,'.dat');
end
CHECK=3;
fullname=fullfile(path,file);
fid = fopen(fullname,'w');
fprintf(fid,'%8.2f %8.2f %8.4f\r\n',XYZ');
fclose(fid);
CHECK=4;


end

