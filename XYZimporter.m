function [GRIDDATA,INFO] = XYZimporter
%XYZimport takes XYZ tabular data and converts it to a matrix usable by
%MATLAB. Also outputs the grid spacing and other pertinent data in the INFO
%struct. Field are xmin, xmax, ymin, ymax, xstep, ystep, xlabel, ylabel.
%Each field is a scalar with the exception of xlabel and ylabel, which are
%vectors.

%pull in data from the file
DATAIN=uiimport;
ff=fields(DATAIN);
FIELDDATA=DATAIN.(ff{1});

%calculate grid extents
xmin=min(FIELDDATA(:,1));
xmax=max(FIELDDATA(:,1));
ymin=min(FIELDDATA(:,2));
ymax=max(FIELDDATA(:,2));

%getting grid cell size
%calculate x spacing
QS=length(FIELDDATA(:,1))-1;
for i = 1:QS
diffx(i)=FIELDDATA(i+1,1)-FIELDDATA(i,1);
end
xstep=min(abs(diffx(abs(diffx)>0)));

%calculate y spacing
QS=length(FIELDDATA(:,2))-1;
for i = 1:QS
diffy(i)=FIELDDATA(i+1,2)-FIELDDATA(i,2);
end
ystep=min(abs(diffy(abs(diffy)>0)));

%calculate grid size
xsize=((xmax-xmin)/xstep)+1;
ysize=((ymax-ymin)/ystep)+1;

%populate blank grid data
GRIDDATA=NaN(ysize,xsize);

%convert XYZ tabular data to array format
for i = 1:length(FIELDDATA(:,1))
    GRIDDATA((FIELDDATA(i,2)-ymin)/ystep+1,(FIELDDATA(i,1)-xmin)/ystep+1)=FIELDDATA(i,3);
end

%flip the data right side up (for plotting in MATLAB 1 starts at the top)
%whereas most geologic softwares plot 1 toward the bottom
GRIDDATA=flip(GRIDDATA);

%create label vectors for plotting in figures
xlabel=xmin:xstep:xmax;
ylabel=ymin:ystep:ymax;

%accumulate grid info
INFO.xmin=xmin;
INFO.xmax=xmax;
INFO.ymin=ymin;
INFO.ymax=ymax;
INFO.xstep=xstep;
INFO.ystep=ystep;
INFO.xlabel=xlabel;
INFO.ylabel=ylabel;

% %Data QC by plotting color figure of data
% clim=[-3 3];
% imagesc(xlabel,ylabel,GRIDDATA,clim);
% pbaspect([1 (ymax-ymin)/(xmax-xmin) 1]);

end

