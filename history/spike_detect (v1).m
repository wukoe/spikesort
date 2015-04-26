% spike detection by threshold 
% DQ= detection quality information {channel}logical[spikeAmt,2]
% FALSE: normal; bit1 TRUE: too close to the previous spike; bit2 TRUE: over-threshold "plateau" too wide; 
% 4: low amplitude?
% 2 modes of output: 
%   spikedetect(X) get the full length [0,1] symbol
%   spikedetect(X,'index') get the index of time of the spikes
%   [SD,Q]=spikedetect(...), Q as quality marks
function [SD,varargout]=spike_detect(X,srate,varargin)
%%% parameter setting - these parameters are ordered according to procssing of data
% Whether to remove the outliner
useOLrm=false;

% Threshold method choice
thresMethod=2;
% parameters related to different choice
if thresMethod==1
    % length of signal to produce moving threshold (s):
    movLen=1; % ! tested
    % moving threshold - above average - by folds of STD
    movThres=4;
elseif thresMethod==2 % Quiroga method
    movLen=10; % length (s) of window 
    movThres=4.2; % folds of estimated STD
elseif thresMethod==3 % RMS calculation
    movLen=2;
    movThres=1.5;
end

% Quality check related
spikeWidthThres=2; % (ms) should be smaller
spikeIntervalThres=3; % (ms) should be bigger


%%%%%%%%%%%%%%%%% data processing
[pntAmt,chAmt]=size(X);
movLen=srate*movLen;

% whether the quality information of detection is needed in output
if nargout==2
    bDInfo=true;
else
    bDInfo=false;
end

% About the output spike train data format
if nargin==1
    bIdxOut=false;
else % index output
    bIdxOut=true;
end

%%% removing low-freq noise & base line of signal & threshold results
% 1) remove outliers - method: detect all points larger than 10 times of
% STD, replace by 0 (this is actually for better STD estimation, STD will be updated later)
X=rmmean(X);
if useOLrm
    ssd=std(X);    
    for chi=1:chAmt
        I=(X(:,chi)>10*ssd(chi)) | (X(:,chi)<-10*ssd(chi));
        X(I,chi)=0;
    end    
end
% % now update the STD
% ssd=std(X);


% 2) detect spikes based on local time window
% Get the windows for local average baseline
[segm,sega]=cutseg(pntAmt,movLen);

%%% Detect all "over-threshold" segments
% the output: array of bites
SD=false(pntAmt,chAmt);
% treat each channel separately
for segi=1:sega
    if thresMethod==3
        % The RMS method need to separate the interval into smaller bins
        [sbseg,sba]=cutseg(segm(segi,2)-segm(segi,1)+1,srate*0.02);
        if sbseg(end,2)-sbseg(end,1)<(sbseg(1,2)-sbseg(1,1))/2
            sbseg(end-1,2)=sbseg(end,2);
            sbseg(end,:)=[];
            sba=sba-1;
        end
    end
    
    for chi=1:chAmt
        binX=X(segm(segi,1):segm(segi,2),chi);
        
        switch thresMethod
            case 1
                % Get the average number without those "over threshold" part
                movAvg=mean(binX);
                ssd=std(binX);
                thresRes=(binX>movAvg+movThres*ssd) | (binX<movAvg-movThres*ssd);

                % Repeat until no point in the pool for calculating moving average
                % is above the threshold
                while multilogic(thresRes,'or')
                    binX(thresRes)=[];
                    if isempty(binX)
                        break
                    end
                    movAvg=mean(binX);
                    ssd=std(binX);
                    thresRes=(binX>movAvg+movThres*ssd) | (binX<movAvg-movThres*ssd);
                end
                
                % restore the "binX" variable
                binX=X(segm(segi,1):segm(segi,2),chi);
                
            case 2
                % Estimation by median
                % see [Quian Quiroga, 2004, Unsupervised Spike Detection and Sorting with Wavelets and Superparamagnetic Clustering]
                movAvg=mean(binX);
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
        
        % Get the overthreshold after processed average - all that matters
        % is the "movAvg" variable
        thresRes=(binX>movAvg+movThres*ssd);        
        SD(segm(segi,1):segm(segi,2),chi)=thresRes;
    end
end

%%%%%%%%%%%%%% mark the exact location of spikes
% Basic strategy: 1. one single spike for each interval; 
% 2. The spike mark is supposed to be at the location of maximum of the interval.
% Currently the spike quality will not affect the detection results,
% rather, it is give in the information data assigned next section.

% save the old data if needed in later quality analysis
if bDInfo
    OSD=SD;
end

% Channel by channel separately
for chi=1:chAmt
    % run through the whole series
    pti=1;
    while pti<=pntAmt
        if SD(pti,chi) % If found the start edge of "plateau" 
            % mark the start edge loc
            ps=pti;
            % search for the end edge
            while pti<=pntAmt && SD(pti,chi)
                pti=pti+1;
            end
            % mark the end edge
            pe=pti-1;
            
            % find the maximum signal loc
            temp=X(ps:pe,chi);
            [~,idx]=max(temp);
            % set all points in this SD interval 0, leaving only the max loc
            SD(ps:pe,chi)=false;
            SD(ps+idx-1,chi)=true;
        end
        pti=pti+1;
    end
end

% if output is set to index mode
if bIdxOut
    I=cell(chAmt,1);
    for chi=1:chAmt
        I{chi}=find(SD(:,chi));
    end
    SD=I;
end

%%%%%%%%%%%%%% output the detection information
% 1. check the width of interval (longer than certain number can't be a real spike) *** not done yet; 
% 2. check the interval between spikes (but this should be considered
% with spike sorting, because 2 neighboring ones may belong to different
% cells !!!)
if bDInfo
    DQ=cell(chAmt,1);
    for chi=1:chAmt
        I=find(SD(:,chi));
        sAmt=length(I);
        DQ{chi}=zeros(sAmt,1,'uint8');
        
        % Width check
        % Transfer unit to (pts)
        spikeWidthThres=round(spikeWidthThres/1000*srate);
        
        tq=false(sAmt,8);
        for si=1:sAmt
            % find the area start and rear ends
            ps=I(si);
            while ps>0 && OSD(ps,chi)
                ps=ps-1;
            end
            ps=ps+1;
            pe=I(si);
            while pe<=pntAmt && OSD(pe,chi)
                pe=pe+1;
            end
            pe=pe-1;
            
            if pe-ps>spikeWidthThres
                tq(si,1)=true;
            end
        end
        
        % Spike interval check - assigned to the rear-end spike
        % Transfer unit to (pts)
        spikeIntervalThres=round(spikeIntervalThres/1000*srate);
        
        for si=2:sAmt
            if I(si)-I(si-1)<spikeIntervalThres
                tq(si,2)=true;
            end
        end
        
        DQ{chi}=logical2num(tq);%%%%%%
    end
    
    varargout{1}=DQ;
end
