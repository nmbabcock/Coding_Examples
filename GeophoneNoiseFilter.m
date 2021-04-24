Function [RAW,RAWC,NOISE,NOISEC,CLEAN,CLEANC] = FilterRT(FF,LF,SL,LT,Vert,Xline,Inline,Mic,Pilot,FTL,ampK,ampFREQ,MULT);


%This program is a simulation of a real-time air noise filter. It uses a
%for loop to pass data point-by-point to the filtering algorithm. This
%filter method modifies the amplitude of a microphone and uses the phase
%from a geophone to create a frequency-domain noise estimate. This is
%inverted back to the time domain to produce a time-domain noise estimate.
%The function is called using
%[RAW,RAWC,NOISE,NOISEC,CLEAN,CLEANC]=
%FilterRT(FF,LF,SL,LT,Vert,Xline,Inline,Mic,Pilot,FTL,ampK,ampFREQ,MULT)
%
%input variables
%FF - First File number
%LF - Last File number
%SL - Sweep length
%LT - listen time
%Vert - Vertical channel number
%Xline - Crossline channel number
%Inline - Inline channel number
%Mic - Microphone channel number
%Pilot - pilot channel number
%FTL - Fourier transform length
%ampK - amplitude filter kernel
%ampFREQ – amplitude kernel frequency axis locations
%MULT - final signal adjustment
%
%output variables
%RAW - uncorrelated data from the seismic record
%RAWC - correlated data from the seismic record
%NOISE – uncorrelated noise estimate for each geophone component
%NOISEC – correlated noise estimate for each geophone component
%CLEAN – uncorrelated noise free geophone record
%CLEANC - correlated noise free geophone record
%
%NOTES
%This program works for a single 3-C geophone and microphone combo. The
%input data (released as Data) can be attached to any channel. The output
%data (NOISE and CLEAN)is organized as such: Vertical-CH1, Crossline-CH2,
%Inline-CH3, Microphone-CH4.
%
%
%The program utilizes the Gabor transform from the Linear
%Time Frequency Analysis Toolbox (Søndergaard et al., 2011)
%[The linear time frequency analysis toolbox: International Journal of Wavelets,
%Multiresolution and Information Processing, 10, website: ltfat.sourceforge.net.],
%and reads SEGY files using SegyMAT (Hansen, 2011).
%[SegyMAT: Revision 1.5, website: segymat.sourcefourge.net.]

%ABOUT THE FILTER
%The filter works by storing the most recent data samples in a buffer. This
%buffer is mirrored and a gaussian window is applied. This action allows
%for the filter to work in real-time. If the buffer was not mirrored the
%filter would work on the central sample, creating a delay of 1/2 the
%Fourier transform length. The microphone data is modified to mimic the
%geophone's amplitude decay (with respect to sound distance). an overview
%of this operation is available within the code. After the symmetric data
%buffer is created a FFT is performed, to convert to the frequency domain.
%In this domain the Microphone is multiplied by each geophones Filter
%Amplitude Kernel. The magnitude of the filtered microphone is used, along
%with the phase from the geophone FFT's to create a frequency-domain noise
%estimate. After conversion back to the time-domain, the central sample of
%the noise estimate (which corresponds to the current time) is subtracted
%from the original data. This produces a single geophone signal-sample
%which should be free of air noise. After this the circular data buffers
%are rotated and the filter process is repeated for the next data sample.


	%pre allocate variables
NOISE=zeros(16000,3);
CLEAN=zeros(16000,3);
MAXr=zeros(LF-FF+1,1);
Mres=zeros(10,1);


TMPA=FF;

for n=1:LF-FF+1 %Start File iteration loop

TMPB=sprintf('%i.sgy',TMPA);
[Data,SegyTraceHeaders,SegyHeader]=ReadSegy(TMPB);
RAW(:,n,:)=Data(:,[Vert Xline Inline Mic Pilot]);
TMPA=TMPA+1;
dt=SegyHeader.time(2)-SegyHeader.time(1);
ts=length(RAW(:,1,1));
xaxisff=(1/dt)*linspace(0,1,2*FTL-1);
Akern=interp1(ampFREQ,ampK,xaxisff,'spline');


	%create circular data buffer for FFT
PastV=zeros(FTL,1);
PastX=zeros(FTL,1);
PastI=zeros(FTL,1);
PastM=zeros(FTL,1);


	%for loop simulates real-time filtering of recorded data
for m=1:length(RAW(:,1,1)); %START MAIN FOR LOOP


	%pass stored data to live variable to mimic real-time filtering.
CdataV=RAW(m,n,1);
CdataX=RAW(m,n,2);
CdataI=RAW(m,n,3);
CdataM=RAW(m,n,4);


	%fill new value in data buffer
PastV(FTL)=CdataV;
PastX(FTL)=CdataX;
PastI(FTL)=CdataI;
PastM(FTL)=CdataM;


%Store local peak amplitudes for distance estimation
%local peak amplitudes are picked for each 1000 samples
%they are stored in a buffer for 10000 samples
%the maximum of the local peaks is used in a distance scaling factor
%this factor adjusts the microphone's 1/r response to sound to mimic
%the 1/r^2 response of the geophone. This is calibrated so that sounds
%coming from 8 meters appear as the same ampltude on both Mic and Geop

M(2)=max(abs(PastM));
if rem(m,1000)==0
 Mres=circshift(Mres,[-1 0]);
 Mres(10)=0;
end

if M(2)>=Mres(10)
 Mres(10)=M(2);
 if M(2)>=max(Mres)
    M(1)=M(2);
  end
end

if M(2)<max(Mres)
 M(1)=max(Mres);
end

if M(2)>MAXr(n)
 MAXr(n)=M(2);
end


	%creation of symmetric signal for real-time FFT
SymV=PastV;
SymX=PastX;
SymI=PastI;


%microphone amplitudes are scaled by M(1)/8.5362e3 which is
%[max(recent peak amplitudes)]/[peak amplitude of sound from 8m]
%This adjusts the microphone to have a psuedo 1/r^2 response
%This is only used on the Vertical geophone component

SymMV=PastM*(M(1)/8.5362e3);

%The Inline and Xline components don't need affected mic inputs
SymM=PastM;


	%fill out symmetric signal
SymV(FTL+1:2*FTL-1)=fliplr(SymV(1:FTL-1));
SymX(FTL+1:2*FTL-1)=fliplr(SymX(1:FTL-1));
SymI(FTL+1:2*FTL-1)=fliplr(SymI(1:FTL-1));
SymM(FTL+1:2*FTL-1)=fliplr(SymM(1:FTL-1));
SymMV(FTL+1:2*FTL-1)=fliplr(SymMV(1:FTL-1));


	%calculate fft of current data buffer
VftS=fft(gausswin(2*FTL-1).*SymV);
XftS=fft(gausswin(2*FTL-1).*SymX);
IftS=fft(gausswin(2*FTL-1).*SymI);
MftS=fft(gausswin(2*FTL-1).*SymM);
MftSV=fft(gausswin(2*FTL-1).*SymMV);


	%calculate phase for each geophone FFT
VphS=angle(VftS);
XphS=angle(XftS);
IphS=angle(IftS);


	%apply amplitude filter kernel to microphone FFT
filtV=(abs(MftSV)).*(Akern(:,1));
filtX=(abs(MftS)).*(Akern(:,2));
filtI=(abs(MftS)).*(Akern(:,3));


	%filter Vertical geophone trace:
%apply geophone phase to filtered microphone amplitude
%to create Noise estimate FFT
RR=abs(filtV).*cos(VphS);
II=abs(filtV).*sin(VphS);
VftNEWS=complex(RR,II);

	%IFFT to create local-noise estimate time signal
TMP4=real(ifft(VftNEWS));
TMP6=TMP4(1:2*FTL-1);
NOISE(m,n,1)=TMP6(FTL);


	%subtract current-time noise estimate from geophone signal
CLEAN(m,n,1)=RAW(m,n,1)-MULT(1)*NOISE(m,n,1);


	%filter Crossline geophone trace
RR=abs(filtX).*cos(XphS);
II=abs(filtX).*sin(XphS);
XftNEWS=complex(RR,II);
TMP4=real(ifft(XftNEWS));
TMP6=TMP4(1:2*FTL-1);
NOISE(m,n,2)=TMP6(FTL);
CLEAN(m,n,2)=RAW(m,n,2)-MULT(2)*NOISE(m,n,2);


	%filter Inline geophone trace
RR=abs(filtI).*cos(IphS);
II=abs(filtI).*sin(IphS);
IftNEWS=complex(RR,II);
TMP4=real(ifft(IftNEWS));
TMP6=TMP4(1:2*FTL-1);
NOISE(m,n,3)=TMP6(FTL);
CLEAN(m,n,3)=RAW(m,n,3)-MULT(3)*NOISE(m,n,3);
if rem(m,1000)==0
	TMP8=sprintf('finished loop #%i',m);
	disp(TMP8)
end

	%advance data buffer 1 time sample
PastV=circshift(PastV,[-1 0]);
PastX=circshift(PastX,[-1 0]);
PastI=circshift(PastI,[-1 0]);
PastM=circshift(PastM,[-1 0]);
end %END MAIN FOR LOOP


	%correlate data
%Raw data
TMP2=xcorr(RAW(:,n,1),RAW(:,n,5),SL*(1/dt));
RAWC(:,n,1)=TMP2(SL*(1/dt)+1:SL*(1/dt)+1+LT*(1/dt));
TMP2=xcorr(RAW(:,n,2),RAW(:,n,5),SL*(1/dt));
RAWC(:,n,2)=TMP2(SL*(1/dt)+1:SL*(1/dt)+1+LT*(1/dt));
TMP2=xcorr(RAW(:,n,3),RAW(:,n,5),SL*(1/dt));
RAWC(:,n,3)=TMP2(SL*(1/dt)+1:SL*(1/dt)+1+LT*(1/dt));
TMP2=xcorr(RAW(:,n,4),RAW(:,n,5),SL*(1/dt));
RAWC(:,n,4)=TMP2(SL*(1/dt)+1:SL*(1/dt)+1+LT*(1/dt));

for o=1:3

	%Noise estimate
	TMP2=xcorr(NOISE(:,n,o),RAW(:,n,5),SL*(1/dt));
	NOISEC(:,n,o)=TMP2(SL*(1/dt)+1:SL*(1/dt)+1+LT*(1/dt));

	%Clean signal
	TMP2=xcorr(CLEAN(:,n,o),RAW(:,n,5),SL*(1/dt));
	CLEANC(:,n,o)=TMP2(SL*(1/dt)+1:SL*(1/dt)+1+LT*(1/dt));
end

TMP=sprintf('Completed file %i',TMPA-1);
disp(TMP)
fclose('all');
end %end file iteration loop