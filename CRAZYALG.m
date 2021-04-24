function [OILIPM, OILIPS, COST, W] = CRAZYALG(OILIP,FOOTAGE)

%This script takes the foot-by-foot steering for each well and the Peak
%Daily Oil IP's to calculate the producability of each zone.

%data is imported in 2 datasets:
%OILIP is a list of daily peak oil IP's by well.
%FOOTAGE is a matrix of the footage drilled in each zone organized by well.

%This algorithm uses a modified solution to solving a system of linear
%equations. It assumes that the IP of a well is a function of the footage
%the well spends in each zone. However, there is a lot of variability, so
%the typical modifiers (e.g. the m and n in Y=mX+nZ) will not be static.
%Instead this method assumes each modifier has a mean and standard
%deviation. The algorithm uses a modified genetic algorithm, with periodic
%mutation, to reach the minimum cost solution.


%pre-define a few indexing variables
m=length(OILIP);
o=length(FOOTAGE(1,:));
n=300;
MNNEW=ones(100,o);
MNLAST=ones(100,o);
MNBACKUP=ones(100,o);
SDVNEW=ones(100,o);
SDVLAST=ones(100,o);
SDVBACKUP=ones(100,o);
COSTLAST=ones(100,1)*999999999999999999999;
idx=1;
MINCOST(1)=999;
TESTPROD=zeros(100,m);

%SET if use Serial of Parallel monte carlo sim
SERIAL=1;

%set number of generations before exit
EX=2000;

%create a matrix of OILIP results to check results against
OILIPM=ones(1,100,m);
for f=1:100
    OILIPM(1,f,:)=OILIP;
end

%the algorithm starts with 100 random distributions
MN=rand(100,o);
SDV=rand(100,o)*.1;

%Run Montecarlo simulation to determine output production of well.
%runs a monte carlo simulation populating the values of the linear
%equations from their associated distributions.
    M=ones(1,o);
    MM=ones(m,o);
    W=ones(100,n,m);
    
startgen=tic;
while min(MINCOST)>100   %LOOP checks error, if low error has been achieved, the loop exits
startloop(idx)=tic;
fprintf('starting generation number %i \n',idx);
newline;

%if SERIAL==1
%This is the serial calculation version of the monte carlo simulation
%STARTSERIAL=tic;
for i=1:n
  
    %iterate for each tested distribution
    for f=1:100
    %for t=1:m  
        %generate string of values from layer distributions
        checkneg=0;
        for j=1:o
            %checkneg makes sure we dont't pull a negative producibility
          while checkneg==0
          tmp=SDV(f,j)*randn+MN(f,j);
            if tmp>=0
              M(j)=tmp;
              checkneg=1;
            end
          end
        checkneg=0;
        end

        %set up matrix form of layer producibility values for dot product
        MM=repmat(M,m,1);

    %calculate dot product for each layer
    %f respresents the tested distribution (1 to 100)
    %i represents the monte carlo sim run number (1 to n)
    %":" is the summed predicted production for the given monte
    %carlo run and test distribution for each well (length m).
    W(f,i,:)=dot(FOOTAGE,MM,2);
   %endre
    end

end
%ENDSERIAL=toc(STARTSERIAL)

%takes the mean and stdev of all monte carlo runs
RESULT(1,:,:)=mean(W,2);
RESULT(2,:,:)=std(W,0,2);
%calculates the fit of each test distribution
%lower values means a closer fit 
COST=RESULT(1,:,:)-OILIPM(1,:,:);

% else
% %this is the parallel version of the monte carlo simulation
% %STARTPAR=tic;
% for f=1:100
% for i=1:m   
%  zum=0;
%     for j=1:o %USE PARFOR???? somehow makes really slow
%         %zum=0;
%        for k=1:n
%             %checkneg makes sure we dont't pull a negative producibility
%             checkneg=0;
%           while checkneg==0
%           tmp=SDV(f,j)*randn+MN(f,j);
%             if tmp>=0
%               zum=zum+tmp*FOOTAGE(i,j);
%               checkneg=1;
%             end
%           end
%        end
%        %TESTPRODZ(f,i,j)=zum/n;
%     end
%   %zum
% TESTPROD(f,i)=zum/n;
%     
% end
% end
% %ENDPAR=toc(STARTPAR)
% %TESTPRODZ
% %TESTPROD=sum(TESTPRODZ,3);
% TMP(1,:,:)=TESTPROD;
% COST=TMP(1,:,:)-OILIPM(1,:,:);
% end


disp('Completed montecarlo simulation')

COSTMAG=ones(1,100);
Q=ones(48,1);
size(Q);
for f=1:100
    size(Q);
    size(COST);
    Q(:)=COST(1,f,:);
COSTMAG(f)=norm(Q);
end
[LOW,I]=sort(COSTMAG,'ascend','ComparisonMethod','abs');

%save lowest error from this round of testing
MINCOST(idx)=LOW(1);
BESTMEAN=MN(I(1),:);
BESTSDV=SDV(I(1),:);

if MINCOST(idx)<100
    break
end

%this genetic algorithm didn't reduce error over time
% %keep the best 5 to preserve genetic material for future runs
% for i=1:5
%     MNNEW(i,:)=MN(I(i),:);
%     SDVNEW(i,:)=SDV(I(i),:);
% end
% i=i+1;
% %Genetic mixing method A:
% %take each parent, swap one of the 19 genes at random with each other
% %parent. Then repeat the process swapping two genes at random.
% 
% %take 10 best test cases and mix them together to create new test cases
%     for j=1:8
%         %loop replaces one gene from parent DNA
%         for k=j:9
%         MNNEW(i,:)=MN(j,:);
%         SDVNEW(i,:)=SDV(j,:);
%         
%         q=randi(19);
%         MNNEW(i,q)=MN(I(k+1),q);
%         SDVNEW(i,q)=SDV(I(k+1),q);
%         i=i+1;
%         end
%         
%         %loop replaces two genes from parent DNA
%         for k=j:9
%         MNNEW(i,:)=MN(j,:);
%         SDVNEW(i,:)=SDV(j,:);
%         
%         q=randi(19);
%         MNNEW(i,q)=MN(I(k+1),q);
%         SDVNEW(i,q)=SDV(I(k+1),q);
%         
%         q=randi(19);
%         MNNEW(i,q)=MN(I(k+1),q);
%         SDVNEW(i,q)=SDV(I(k+1),q);
%         i=i+1;
%         end       
%     end
% 
% %add some new test datasets to the pool, based off the average gene values
% %and standard deviations
% for i=94:100
%     for k=1:19
%        checkneg=0;
%        while checkneg==0
%        tmp=std(MN(k,:))*randn+mean(MN(k,:));
%           if tmp>=0
%           MNNEW(i,k)=tmp;
%           checkneg=1;
%           end
%        end
%        
%        checkneg=0;
%        while checkneg==0
%        tmp=std(SDV(k,:))*randn+mean(SDV(k,:));
%           if tmp>=0
%           SDV(i,k)=tmp;
%           checkneg=1;
%           end
%        end
%     end
% end



%pull the MN and SDV from backup, to be used in genetic algorithm
MNLAST=MNBACKUP;
SDVLAST=SDVBACKUP;
%save current MN and SDV to backup to be used next iteration
MNBACKUP=MN;
SDVBACKUP=SDV;

%Genetic mixing method B:
%If a child performed better than before, its DNA gets slightly perturbed
%if it did not, the parent gets perturbed again

for i=1:100
    if COSTMAG(i)<COSTLAST(i)
        %disp('new magnitude is better!');
        tmp=randi(19);
        flip=randi(2);
    %change either the mean or standard deviation of one element within the parant    
      if flip==1
        checkneg=0;
        while checkneg==0
            mean(MN(i,:));
         %tmpb=MN(i,tmp)+mean(MN(i,:))/100*randn;
         tmpb=MN(i,tmp)+0.1*randn;
         if tmpb>=0
          MN(i,tmp)=tmpb;
          checkneg=1;
         end
        end
      else
        checkneg=0;
        while checkneg==0
         %tmpb=SDV(i,tmp)+mean(SDV(i,:))/100*randn;
         tmpb=SDV(i,tmp)+0.1*randn;
          if tmpb>=0
           SDV(i,tmp)=tmpb;
           checkneg=1;
          end
         end
      end
    else
        tmp=randi(19);
        flip=randi(2);
    %change either the mean or standard deviation of one element within the parant    
      if flip==1
        checkneg=0;
        while checkneg==0
         %tmpb=MNLAST(i,tmp)+mean(MNLAST(i,:))/100*randn;
         tmpb=MNLAST(i,tmp)+0.1*randn;
         if tmpb>=0
          MN(i,tmp)=tmpb;
          checkneg=1;
         end
        end
      else
        checkneg=0;
        while checkneg==0
         %tmpb=SDVLAST(i,tmp)+mean(SDVLAST(i,:))/100*randn;
         tmpb=SDVLAST(i,tmp)+0.1*randn;
          if tmpb>=0
           SDV(i,tmp)=tmpb;
           checkneg=1;
          end
         end
      end
    end
end
COSTLAST=COSTMAG;


%MN=MNNEW;
%SDV=SDVNEW;
fprintf('current minimum error is %d barrels \n',MINCOST(idx))
newline;
endloop(idx)=toc(startloop(idx));
fprintf('loop duration: %d \n',endloop(idx))
newline;
idx=idx+1;

if idx==EX+1
    MINCOST(idx)=0;
end
end
endgen=toc(startgen);

fprintf('The algorithm has completed after running through %i generations in %d seconds \n',idx,endgen)
newline;
fprintf('The minimum achieved arror is %d barrels \n',min(MINCOST))
newline;
OILIPM=MN(I(1),:);
OILIPS=SDV(I(1),:);
COST=MINCOST;

%to do:
%---DONE---   create montecarlo simulator
%---DONE---   calculate cost of function results (results-OILIP => 0)
%---DONE---   choose best results
%---DONE---   merge best results and create child test cases
%---DONE---   add in mutation clause
%---DONE---   track minimum error per generation
%---DONE---   create abandonment case



end



  




