% spike detection algorithm
%   [SD,SQ,SA]=spikedetect(X,srate) 
% output the spike data SD (logical format{0,1} symbol), quality and
% amplitude of spikes (SA).
%   [SD,SQ,SA,rawSQ]=spikedetect(...), give raw SQ before filtering.
% SQ= detection quality information {channel}logical[spikeAmt,2]
% FALSE: normal; bit1 TRUE: too close to the previous spike; bit2 TRUE: over-threshold "plateau" too wide; 
% 4: low amplitude?
function [SD,SQ,SA,SW,varargout]=spike_detect(X,srate,paras)
%%%%%%%%%%%%%% Option parameter default setting (ordered according to procssing of data)
%%% Default
% Whether to remove the outliner
useOLrm=false;
OLrmThres=150;% threshold for removal: fold of STD in a bin

% Threshold method choice
thresMethod=2;
% parameters related to different choice. MovLen: the segment for
% calculating the local average.
if thresMethod==1 % Basic method
    % length of signal to produce moving threshold (s):
    movLen=10; % ! tested
    % moving threshold - above average - by folds of STD
    movThresPos=4;
elseif thresMethod==2 % Quiroga method
    movLen=1; % length (s) of window   <<<<<<<<<<<<%%%%%%%%%%%%%
    movThresPos=5; % folds of estimated STD % 5 also works I think
    movThresNeg=5; % 
elseif thresMethod==3 % RMS calculation
    movLen=2;
    movThresPos=1.5;
end

% ɸѡspike
% �Ƿ��Ƴ������ڷ�Χ֮���spike
bBadWidthRemove=true;
spikeWidthThres=[0,0.8]; % (ms) should be in between
% * ����ע�⣬ʵ�����ֺܶ��׼��spike����ֻ��1�������㣨��Ӧ<0.1ms��

% �Ƿ�ϲ�һ��complex����Ķ����
bMergePeak=true;
ppThres=1; % maximum Peak-peak distance allowed. *����Ҫ��refractory period ��һ�£���ΪĿǰ����������spike������Դ�ڲ�ͬ��ϸ����
% totalLenThres=4; % maximum total length threshold.
% PNpeakRatioThres=0.5;

% �Ƿ�������������ȷ��λ�������λ��
bAccuPeak=false;
apSearchIntH=0.2; %(ms) -half of the interval to do spline (extend left and right from the maximum location)

% SQ�����¼interval��С����ֵ
spikeIntervalThres=5; % (ms) should be bigger than
% * spikeIntervalThres Ӧ�ñ�postww or preww������������smoothing other spike ���ṩ����
% ��������DQ��spikeIntervalThres��Ŀ�м�¼Ϊ����������

%%% User input
if ~isempty(paras)
    if isfield(paras,'movThres')
        movThresPos=paras.movThres;
        movThresNeg=movThresPos;
    elseif isfield(paras,'')
        
    end
end
% Keep or disregard Positive/Negative peaks
bNP=paras.bNP;
bPP=paras.bPP;


%%% Process
pntAmt=length(X);
movLen=srate*movLen;

% Transfer unit in (ms) to (pts)
spikeWidthThres=round(spikeWidthThres/1000*srate);
spikeIntervalThres=round(spikeIntervalThres/1000*srate);
ppThres=round(ppThres/1000*srate);
% totalLenThres=round(totalLenThres/1000*srate);
apSearchIntH=round(apSearchIntH/1000*srate);


%%%%%%%%%%%%%%% Mark all "over-threshold" segments
% the output: array of bites
NSD=false(pntAmt,1);
PSD=NSD;

% Get the windows for bin - detect spikes based on local time window
[segm,sega]=cutseg([1,pntAmt],movLen);
for segi=1:sega
    % The RMS method need to separate the interval into smaller bins (20ms)
    % - determine the small bin segmentation.
    if thresMethod==3        
        [sbseg,sba]=cutseg([1,segm(segi,2)-segm(segi,1)+1],srate*0.02);
        if sbseg(end,2)-sbseg(end,1)<(sbseg(1,2)-sbseg(1,1))/2
            sbseg(end-1,2)=sbseg(end,2);
            sbseg(end,:)=[];
            sba=sba-1;
        end
    end
    
    
    binX=X(segm(segi,1):segm(segi,2));

    % Remove outliers - method: detect all points larger than 15 times (OLrmThres) of
    % STD, replace by 0 (STD might need to be updated later)
    if useOLrm
        binX=rmmean(binX);
        ssd=std(binX);           
        I=(binX>OLrmThres*ssd) | (binX<-OLrmThres*ssd);
        %<<< I need to be extended
        binX(I)=0;            
    end

    %%% Determine the mean and ssd estimation for each method.
    switch thresMethod
        case 1
            % Get the average number without those "over threshold" part
            movAvg=mean(binX);
            ssd=std(binX);
            thresRes=(binX>movAvg+movThresPos*ssd) | (binX<movAvg-movThresPos*ssd);

            % Repeat until no point in the pool for calculating moving average
            % is above the threshold
            while multilogic(thresRes,'or')
                binX(thresRes)=[];
                if isempty(binX)
                    break
                end
                movAvg=mean(binX);
                ssd=std(binX);
                thresRes=(binX>movAvg+movThresPos*ssd) | (binX<movAvg-movThresPos*ssd);
            end

            % restore the "binX" variable
            binX=X(segm(segi,1):segm(segi,2));

        case 2
            % Estimation by median
            % see [Quian Quiroga, 2004, Unsupervised Spike Detection and Sorting with Wavelets and Superparamagnetic Clustering]
            movAvg=mean(binX);
            binX=binX-movAvg; %%%%
            ssd=median(abs(binX))/0.6745;

        case 3
            binNum=zeros(sba,1);
            for sbi=1:sba
                binNum(sbi)=mean(binX(sbseg(sbi,1):sbseg(sbi,2)).^2); % mean square
            end
            binNum=sort(binNum,'ascend');
            ssd=mean(binNum(6:25));
            movAvg=mean(binX);
    end

    %%% Get the over-threshold part after processed average - all that matters
    % is the "movAvg" variable.
    % Positive peaks.
    if bPP
        thresRes=(binX>movThresPos*ssd);        
        PSD(segm(segi,1):segm(segi,2))=thresRes;
    end
    % Negative peaks.
    if bNP
        thresRes=(binX<-movThresNeg*ssd);        
        NSD(segm(segi,1):segm(segi,2))=thresRes;
    end
end


%%%%%%%%%%%%%% mark the exact location of spikes
% Basic strategy: 1. one single spike for each interval; 
% 2. The spike mark is supposed to be at the location of maximum of the interval.
% Currently if refractory is set, than following spikes within interval threshold will be ignored.

% Initiate The Spike data quality information: one cell for each channel.
% * It include 3 columns: 
% 1. marke polarity of spike: positive(1)/negative(0); 
% 2. (reserved)
% 3&4. check peak width (larger than certain number can't be a real spike): too wide(1)/normal(0);  
% 5. check the interval between spikes: too close to previous(1)/normal(0)
% (but this later should be considered with spike sorting, because 2 neighboring ones may belong to different cells)
sqBit=5; % total number of SQ flags.
sqPolLoc=1;
sqWidNLoc=3;
sqWidWLoc=4;
sqIntLoc=5;

SD=zeros(0,1);
SA=zeros(0,1); SW=SA;
SQ=false(0,sqBit);
sPA=0;

%%% The positive peaks
if bPP
P=continuous_segment(PSD);
sPA=size(P,1);

% Init SD, SA and SQ
SD=zeros(sPA,1);
SA=zeros(sPA,1); SW=SA;
SQ=false(sPA,sqBit);
% Polarity information in SQ
SQ(:,sqPolLoc)=true;

for si=1:sPA
    % Find the maximum signal location
    p=P(si,:);
    temp=X(p(1):p(2));
    [SA(si),idx]=max(temp);
    SD(si)=p(1)+idx-1;
end
% Peak width check
SW=diff(P,[],2);
SQ(:,sqWidNLoc)=(SW<spikeWidthThres(1));
SQ(:,sqWidWLoc)=(SW>spikeWidthThres(2));
end

%%% The negative peaks
if bNP
P=continuous_segment(NSD);
sNA=size(P,1);

% extend the SD and SQ
SD=[SD; zeros(sNA,1)];
SA=[SA; zeros(sNA,1)];
SQ=[SQ; false(sNA,sqBit)];
% Polarity information in SQ
SQ(sPA+1:sPA+sNA,sqPolLoc)=false;

for si=1:sNA   
    % Find the negative maximum signal location
    p=P(si,:);
    temp=X(p(1):p(2));
    [SA(sPA+si),idx]=min(temp);
    % set all points in this SD interval 0, leaving only the max loc
    SD(sPA+si)=p(1)+idx-1;
end
% Peak width check
tp=diff(P,[],2);
SQ(sPA+1:sPA+sNA,sqWidNLoc)=(tp<spikeWidthThres(1));
SQ(sPA+1:sPA+sNA,sqWidWLoc)=(tp>spikeWidthThres(2));
SW=[SW;tp];
end

%%% Sort spikes list position by time
[SD,I]=sort(SD,'ascend');
SA=SA(I); SW=SW(I);
SQ=SQ(I,:);

%%% To output ORIGINAL SQ.
if nargout==4
    varargout{1}=SQ;
end
clear PSD NSD


%%%%%%%%%%%%%%%%%%%% Checking - remove potentially false-positive spikes
%%% �Ƴ������ڷ�Χ֮���spike
if bBadWidthRemove
    rmlist=find(SQ(:,sqWidNLoc));        
    SD(rmlist)=[];
    SA(rmlist)=[]; SW(rmlist)=[];
    SQ(rmlist,:)=[];
    rmlist=find(SQ(:,sqWidWLoc));        
    SD(rmlist)=[];
    SA(rmlist)=[]; SW(rmlist)=[];
    SQ(rmlist,:)=[];
end

%%% spike���ж�����ʱ�򣬽��кϲ�
if bMergePeak
    % Find continuous series of peaks (peak complex) whose interval < threshold, until
    % reaches maximum cumulation threshold
    D=diff(SD);
    cplxList=continuous_segment(D<=ppThres);
    cplxList(:,2)=cplxList(:,2)+1;

    %%% In each complex, determine the peak with maximum amplitude, and remove the others.
    onetwoRatio=2;
    cplxAmt=size(cplxList,1);
    if cplxAmt>0
        rmlist=[];
        for cpi=1:cplxAmt
            % Find amplitude of peaks in complex
            temp=SA(cplxList(cpi,1):cplxList(cpi,2));
            cplxpn=length(temp);            
            if cplxpn>3
                [atp,idx]=sort(abs(temp),'descend');
                % Ӧ�ù���1�����������Ķ����ʱ����������ķ�ֵ�ر�󣬷�����������complex
                if atp(1)>=atp(2)*onetwoRatio;
                    seleidx=idx(1);
                else
                    seleidx=[]; % ��ѡ��complex�е��κ�
                end
            elseif cplxpn==3
                % Find whether rule2 applies
                ps=sign(temp); % polarity sign
                if isequal(ps,[-1;1;-1]) || isequal(ps,[1;-1;1]) % ����Ƿ�NPN or PNP
                    % Ӧ�ù���2������(PNP or NPN�������м�ķ�Ϊ׼��
                    seleidx=2;
                else
                    % 1/N
                    [~,seleidx]=max(abs(temp));
                    
                    % 2/N
%                     I=find(temp<0,1);
%                     if ~isempty(I)
%                         [~,seleidx]=min(temp);
%                     else
%                         [~,seleidx]=max(temp);
%                     end
                    
                    % 3/N
%                     sps=sum(ps==-1); % number of negative peaks
%                     if sps==0
%                         % ���ж�Ϊ��,�Ե�һ�����ֵ�Ϊ׼
%                         seleidx=1;
%                     elseif sps==1 || sps==2
%                         % Ӧ�ù���3�������ڸ��壬���ǵ�һ���������Ҵ��ڵڶ��������������Ը�Ϊ׼��
%                         [atp,idx]=sort(abs(temp),'descend');
%                         if ps(idx(1))==1 && atp(1)>=atp(2)*onetwoRatio
%                             seleidx=idx(1);
%                         else                            
%                             I=(ps==-1); idx=idx(I);
%                             seleidx=idx(1);
%                         end
%                     else % sps==3 Ӧ�ù���4�����ж���ĸ���
%                         [~,seleidx]=min(temp);
%                     end
                    %e/N
                end                    
            else  % cplxpn==2
                [~,seleidx]=max(abs(temp));
            end
            
            % Remove all others in the complex except the selected
            tp=(cplxList(cpi,1):cplxList(cpi,2))';
            tp(seleidx)=[];
            rmlist=[rmlist;tp];
        end
        
        % Delete the others
        SD(rmlist)=[];
        SA(rmlist)=[]; SW(rmlist)=[];
        SQ(rmlist,:)=[];
    end
end

%%% Apply accurate localization of peak  - by spline modeling  <<< ʩ����
% * ���ﲻ��Ҫ����������ܽӽ��������ͨ������searchInt�ķ�����֤��
if bAccuPeak
    ASD=SD;    
    % 1/N
    xpos=-apSearchIntH:apSearchIntH;
    xx=linspace(-apSearchIntH,apSearchIntH, apSearchIntH*10);
    for k=1:length(SD)
        binX=X(SD(k)-apSearchIntH : SD(k)+apSearchIntH);        
        bl=myspline(xpos,binX,xx);
        if SQ(k,sqPolLoc)
            [~,idx]=max(bl);
        else
            [~,idx]=min(bl);
        end
        % assign new location
        ASD(k)=SD(k)-apSearchIntH-1+xx(idx);
    end
%     % 2/N
%     for k=1:length(SI)         
%             binX=X(SI(k)-apSearchIntH : SI(k)+apSearchIntH);
%             binX=diff(binX);
%             [~,idx]=max(abs(binX));
%             SD(k)=SI(k)-apSearchIntH-1+idx;       
%     end
    % e/N
end

%%% To write the interval information in SQ
D=diff(SD);
I=(D<spikeIntervalThres);
SQ(I,sqIntLoc)=true;
        
end