function [GRIDOUT] = RMS2D(INPUTDATA,MOD,RANGE)
%Calculates fracture intensity from minimum curvature. Could also be used
%to calculate aggregate area values for spacially varying datasets, like
%area weighted porosity around a frac stage, etc.
%   Input arguments are 
%   INPUTDATA: data to be worked on in matrix form.
%   MOD: weighting function, G for gaussian, N for none, may add others.
%   SHAPE: What type of shape to use, S for square, C for circle, others
%   Range: 2 element vector (X,Y) containing the offset in grid elements
%   it is best to convert from distance to grid bins before input.

WORKDATA=nan(size(INPUTDATA)+fliplr(RANGE));
for i  = RANGE(2)+1 : length(WORKDATA(:,1))-RANGE(2)-1
   for j = RANGE(1)+1 : length(WORKDATA(1,:))-RANGE(1)-1
       WORKDATA(i,j)=INPUTDATA(i-RANGE(2),j-RANGE(1));
   end
end
max(max(WORKDATA))
GRIDOUT=nan(size(INPUTDATA));
%idx=1; %QC
%imagesc(INPUTDATA,[-3 3]);idx=idx+1; %QC
%figure(idx);imagesc(WORKDATA,[-3,3]); idx=idx+1; %QC

%set-up modifier kernel
if MOD == 'G'
 %2D Gausian Kernel   
    TMPSTEPA=3/RANGE(1);
    TMPSTEPB=3/RANGE(2);
   for i = 1:2*RANGE(1)+1
       TMPA(i,:)=normpdf([-3:TMPSTEPA:3]);
   end
   
   for i = 1:2*RANGE(2)+1
       TMPB(:,i)=normpdf([-3:TMPSTEPB:3]);
   end
   
   DATAMOD=TMPA.*TMPB;
   DATAMOD=DATAMOD./(max(max(DATAMOD)));
else
    %default case, no modifier
    DATAMOD=ones(fliplr(RANGE)); 
end


%figure(idx);imagesc(DATAMOD);idx=idx+1; %QC
%square the data, S part of RMS
WORKDATA=WORKDATA.*WORKDATA;
%figure(idx);imagesc(WORKDATA);idx=idx+1; %QC
%calculate full RMS for each data point
for i  = RANGE(2)+1 : length(WORKDATA(:,1))-RANGE(2)-1
    i
   for j = RANGE(1)+1 : length(WORKDATA(1,:))-RANGE(1)-1
       if isnan(WORKDATA(i,j)) == 0
            GRIDOUT(i-RANGE(2),j-RANGE(1))=sqrt(nanmean((WORKDATA(i-RANGE(2):i+RANGE(2),j-RANGE(1):j+RANGE(1)).*DATAMOD),'all'));
       else
           GRIDOUT(i-RANGE(2),j-RANGE(1))=NaN(1);
       end
   end
end
%figure(idx);imagesc(GRIDOUT);idx=idx+1; %QC
end

