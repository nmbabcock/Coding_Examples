function [FFTOUT,VECTORTHETATRIMRAD] = FRACAZ(SIZE)
%Calculates Fracture azimuth rose diagrams from curvature

[GRIDDATA,INFO] = XYZimporter;
%accumulate grid info
%INFO.xmin=xmin;
%INFO.xmax=xmax;
%INFO.ymin=ymin;
%INFO.ymax=ymax;
%INFO.xstep=xstep;
%INFO.ystep=ystep;
%INFO.xlabel=xlabel;
%INFO.ylabel=ylabel;

%convert SIZE from feet units to cell units
%Note: X->Columns, Y->Rows
SIZEX=floor(SIZE/INFO.xstep)
SIZEY=floor(SIZE/INFO.ystep)


[ROWS,COLUMNS]=size(GRIDDATA)
%convert any NaN to 0
GRIDDATA(isnan(GRIDDATA))=0;
%pad data by SIZE
ROWSPAD=ROWS+2*SIZEY+(SIZEY-mod(ROWS,SIZEY))
COLUMNSPAD=COLUMNS+2*SIZEX+(SIZEX-mod(COLUMNS,SIZEX))
%calculate the number of FFT calculation points.
%NOTE: this method does not have overlap between evaluation areas
ROWSITERATE=ROWSPAD/SIZEY
COLUMNSITERATE=COLUMNSPAD/SIZEX
WORKDATA=zeros(ROWSPAD,COLUMNSPAD);
for i=1:ROWS
    for j=1:COLUMNS
        WORKDATA(i+SIZEX+floor(mod(ROWS,SIZEX)/2),j+SIZEY+floor(mod(COLUMNS,SIZEY)/2))=GRIDDATA(i,j);
    end
end
size(GRIDDATA)
size(WORKDATA)

if mod(SIZEY,2)==1
FFTOUT=zeros(ROWSITERATE-1,COLUMNSITERATE-1,SIZEY,SIZEX);
else
FFTOUT=zeros(ROWSITERATE-1,COLUMNSITERATE-1,SIZEY+1,SIZEX+1);
end

%Calculate sectored FFTs. This method does not have overlaps
for i=1:ROWSITERATE-1
    for j=1:COLUMNSITERATE-1
        TMP=fft2(WORKDATA(i*SIZEY-floor(SIZEY/2):i*SIZEY+floor(SIZEY/2),j*SIZEX-floor(SIZEX/2):j*SIZEX+floor(SIZEX/2)));
        FFTOUT(i,j,:,:)=TMP;
    end
end

%Calculate Azimuth (THETA) and Offset (RHO) attribute
TMP(:,:)=FFTOUT(1,1,:,:);
GRIDSIZE=size(TMP)
if mod(GRIDSIZE(:,1),2)==0
    TOP=GRIDSIZE(:,1)/2;
    BOTTOM=TOP-1;
else
    TOP=(GRIDSIZE(:,1)-1)/2;
    BOTTOM=TOP;
end

if mod(GRIDSIZE(1,:),2)==0
    LEFT=GRIDSIZE(1,:)/2;
    RIGHT=LEFT-1;
else
    LEFT=(GRIDSIZE(1,:)-1)/2;
    RIGHT=LEFT;
end

for i=-TOP:BOTTOM %number of columns
    for j=-LEFT:RIGHT %number of rows
       [theta,rho] = cart2pol(i,j);
       GRIDRHO(i+TOP+1,j+LEFT+1)=rho;
       GRIDTHETA(i+TOP+1,j+LEFT+1)=theta*(180/pi);
    end
end
size(GRIDTHETA)
size(GRIDRHO)
VECTORTHETA=reshape(GRIDTHETA,[],1);
VECTORRHO=reshape(GRIDRHO,[],1);

edges=0:0.0698131701:2*pi;
centers=2:4;358;
for i=1:ROWSITERATE-1
    for j=1:COLUMNSITERATE-1
        TMPB(:,:)=FFTOUT(i,j,:,:);
        TMPC=fftshift(TMPB);
        FFTOUTSHIFT(i,j,:,:)=TMPC;
        TMPD=reshape(abs(TMPC),[],1);
        VECTORFFTOUT(i,j,:)=TMPD;
        
        VECTORTHETATRIMRAD(i,j,:)=(VECTORTHETA+180)*(pi/180);
        for k=1:length(TMPD)
            if VECTORFFTOUT(i,j,k)<0.01
            VECTORTHETATRIMRAD(i,j,k)=NaN;
            end
        end
        
        VECTORTHETATRIMRADBINNED(i,j,:)=histcounts(VECTORTHETATRIMRAD(i,j,:),edges);
        
    end
end

RSCALE=max(max(max(VECTORTHETATRIMRADBINNED)))

size(VECTORTHETA)
size(VECTORRHO)
size(VECTORFFTOUT)
%TMPE(:)=VECTORFFTOUT(13,7,:);
%figure(1);plot3(VECTORTHETA,VECTORRHO,TMPE,'.');

%l=1;
%figure(2)
%for i=1:13
%    for j=1:25
%        TMPF(:)=VECTORTHETATRIMRAD(i,j,:);
%        subplot(13,25,l),rose(TMPF,90);
%        l=l+1;
%    end
%end

H=(ROWSITERATE-1);
W=(COLUMNSITERATE-1);


%{
set(gcf, 'PaperPosition', [0 0 H W])    % can be bigger than screen 

clims = [-.0003 0.0003];
imagesc(GRIDDATA,clims);
hold;
TMPF(:)=VECTORTHETATRIMRAD(13,7,:);
subplot(H,W,1);polarhistogram(VECTORTHETATRIMRAD(13,7,:),edges);
ax=gca;
ax.RLim = [0 RSCALE];
ax.RTick = [];
ax.ThetaTick = [];
ax.RColor = 'white';
ax.ThetaColor = 'white';
ax.ThetaZeroLocation = 'top';
%polaraxes('ThetaZeroLocation','top','ThetaColor','none','RColor','none')


print(gcf, 'MyFigure.png', '-dpng', '-r300' );   %save file as PNG w/ 300dpi
%}



%convert pixels to inches for axes
xlabel=0:(COLUMNS/100)/COLUMNS:COLUMNS/100;
ylabel=0:(ROWS/100)/ROWS:ROWS/100;

set(gcf,'visible','off','units','inches','PaperPosition', [0 0 COLUMNS/100 ROWS/100])    % can be bigger than screen 
clims = [-.0003 0.0003];
h = axes('units','inches','Position',[0 0 COLUMNS/100 ROWS/100]);
%axes('units','inches','position',[0 0 COLUMNS/100 ROWS/100]);
imagesc(h,GRIDDATA,clims);

%{
TMPF(:)=VECTORTHETATRIMRAD(13,7,:);
subplot(H,W,1);polarhistogram(VECTORTHETATRIMRAD(13,7,:),edges);
ax=gca;
ax.RLim = [0 RSCALE];
ax.RTick = [];
ax.ThetaTick = [];
ax.RColor = 'white';
ax.ThetaColor = 'white';
ax.ThetaZeroLocation = 'top';
%polaraxes('ThetaZeroLocation','top','ThetaColor','none','RColor','none')
%}

%'-r300'

print(gcf, 'MyFigure.png', '-dpng', '-r300');   %save file as PNG w/ 300dpi


