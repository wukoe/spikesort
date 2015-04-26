% spike detection algorithm
%   [SD,A]=spikedetect(X) 
% output the spike data SD (logical format{0,1} symbol), and windowed spike
% signals (for the later feature extraction)
%   [SD,A,DQ]=spikedetect(...), Q as quality marks
% DQ= detection quality information {channel}logical[spikeAmt,2]
% FALSE: normal; bit1 TRUE: too close to the previous spike; bit2 TRUE: over-threshold "plateau" too wide; 
% 4: low amplitude?
function [SD,A,varargout]=spike_detect(X,srate,varargin)
%%% Option parameter default setting - ordered according to procssing of data
% Whether to remove the outliner
useOLrm=true;
OLrmThres=20;% threshold for removal

% Whether to detect negative peaks.
bNegPeak=true;

% Threshold method choice
thresMethod=2;
% parameters related to different choice. MovLen: the segment for
% calculating the local average.
if thresMethod==1 % Basic method
    % length of signal to produce moving threshold (s):
    movLen=1; % ! tested
    % moving threshold - above average - by folds of STD
    movThresPos=4;
elseif thresMethod==2 % Quiroga method
    movLen=1; % 10 % length (s) of window   %%%%%%%%%%%%%%%%%%%%%%%%%
    movThresPos=4.2; % folds of estimated STD
    movThresNeg=-4.2;
elseif thresMethod==3 % RMS calculation
    movLen=2;
    movThresPos=1.5;
end

% 是否移除跟随有巨大负峰的spike
bCheckNegaFollow=false;
checkNegaInt=1.5; % (ms) - within this window after peak

% ***************
bCheckStatistic=false;

% 是否加上不应期（refractory period）
useRefra=true; % The refractory time is the same as spikeIntervalThres
refraInterval=1.5; % (ms) 

% 是否用样条方法精确定位尖峰所在位置
bAccuPeak=true;
apSearchIntH=0.5; % half of the interval to search (left and right)

% window for spike waveform
preww=2.5; % pre-event window width (ms)
postww=5; % post-event

% 是否平滑主峰两侧的其他spike造成的峰
bSmoothOthers=true;

% Quality check related
spikeWidthThres=2; % (ms) should be smaller than
spikeIntervalThres=6; % (ms) should be bigger than
% * spikeIntervalThres 应该比postww or preww大，这样才能在smoothing other spike 中提供保护
% （后者以DQ中spikeIntervalThres项目有记录为激活条件）

%%% Basic information
[pntAmt,chAmt]=size(X);
movLen=srate*movLen;

% If the quality information of detection is needed in output, or need to
% smooth the peaks in a window
if nargout==3 || bSmoothOthers
    nvarargout=nargout-2;
    bDQInfo=true;
else
    bDQInfo=false;
end

% Further processing option parameters
% Transfer unit in (ms) to (pts)
spikeWidthThres=round(spikeWidthThres/1000*srate);
spikeIntervalThres=round(spikeIntervalThres/1000*srate);
refraInterval=round(refraInterval/1000*srate);
checkNegaInt=round(checkNegaInt/1000*srate);
apSearchIntH=round(apSearchIntH/1000*srate);

preww=fix(preww/1000*srate);
postww=fix(postww/1000*srate);


%%%%%%%%%%%%%%% Detect all "over-threshold" segments
% the output: array of bites
SD=false(pntAmt,chAmt);
if bNegPeak
    NSD=SD;
end

% Get the windows for bin - detect spikes based on local time window
[segm,sega]=cutseg(pntAmt,movLen);

movAvgSave=zeros(sega,chAmt);
for segi=1:sega
    % The RMS method need to separate the interval into smaller bins (20ms)
    if thresMethod==3        
        [sbseg,sba]=cutseg(segm(segi,2)-segm(segi,1)+1,srate*0.02);
        if sbseg(end,2)-sbseg(end,1)<(sbseg(1,2)-sbseg(1,1))/2
            sbseg(end-1,2)=sbseg(end,2);
            sbseg(end,:)=[];
            sba=sba-1;
        end
    end
    
    for chi=1:chAmt
        binX=X(segm(segi,1):segm(segi,2),chi);
        
        % Remove outliers - method: detect all points larger than 15 times (OLrmThres) of
        % STD, replace by 0 (STD might need to be updated later)
        if useOLrm
            binX=rmmean(binX);
            ssd=std(binX);           
            I=(binX>OLrmThres*ssd) | (binX<-OLrmThres*ssd);
            binX(I)=0;            
        end
        
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
        
        % Get the over-threshold part after processed average - all that matters
        % is the "movAvg" variable
        thresRes=(binX>movAvg+movThresPos*ssd);        
        SD(segm(segi,1):segm(segi,2),chi)=thresRes;
        if bNegPeak
            thresRes=(binX<movAvg+movThresNeg*ssd);        
            NSD(segm(segi,1):segm(segi,2),chi)=thresRes;
        end
        
        % Save the movAvg, could be needed (like negative detection)
        movAvgSave(segi,chi)=movAvg;
    end
end


%%%%%%%%%%%%%% mark the exact location of spikes
% Basic strategy: 1. one single spike for each interval; 
% 2. The spike mark is supposed to be at the location of maximum of the interval.
% Currently if refractory is set, than following spikes within interval threshold will be ignored.

% save the old data if needed in later quality analysis
if bDQInfo
    OSD=SD;
    if bNegPeak
        ONSD=NSD;
    end
end

%%% if the postive/negative issue need to be considered, add about here %%%%%%%%%%%%%%

for chi=1:chAmt
    %%% Run through the whole series
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
            
            % find the maximum signal location
            temp=X(ps:pe,chi);
            [~,idx]=max(temp);
            % set all points in this SD interval 0, leaving only the max loc
            SD(ps:pe,chi)=false;
            SD(ps+idx-1,chi)=true;
        end
        pti=pti+1;
    end
    % Negative
    if bNegPeak
    while pti<=pntAmt
        if NSD(pti,chi) % If found the start edge of "plateau" 
            % mark the start edge loc
            ps=pti;
            % search for the end edge
            while pti<=pntAmt && NSD(pti,chi)
                pti=pti+1;
            end
            % mark the end edge
            pe=pti-1;
            
            % find the maximum signal location
            temp=X(ps:pe,chi);
            [~,idx]=min(temp);
            % set all points in this SD interval 0, leaving only the max loc
            NSD(ps:pe,chi)=false;
            NSD(ps+idx-1,chi)=true;
        end
        pti=pti+1;
    end
    end
    
    %%% Apply negative following check 
    if bCheckNegaFollow
        SI=find(SD(:,chi));
        for k=1:length(SI)
            si=SI(k);
            % determine the segment of this point - in order to get the seg
            % baseline stored ahead.
%             Ilow=(si>=segm(:,1));
%             Ihigh=(si<=segm(:,2));
%             I=Ilow & Ihigh;
%             idx=find(I);            
%             bl=movAvgSave(idx,chi);

            rs=round(5/1000*srate);
            temp=X(si-rs:si+2*rs,chi);
            bl=mean(temp);
            
            % min value in the following period.
            [~,tp]=min(X(si:si+checkNegaInt,chi));
            if bl-X(si+tp,chi) > 2*(X(si,chi)-bl)
                SD(si,chi)=false;
            end
        end
        
        % negative peak here **************
    end
    
    %%% Apply statistical check
    if bCheckStatistic
        %
    end
    
    %%% Apply refractory period
    if useRefra 
        SI=find(SD(:,chi));
        SIdif=diff(SI); % SIdif = interval
        I=(SIdif<refraInterval); % find interval below spikeIntervalThres
        SD(SI(I+1),chi)=false; % set coresponding to 0
        if bNegPeak
            SI=find(NSD(:,chi));
            SIdif=diff(SI); % SIdif = interval
            I=(SIdif<refraInterval); % find interval below spikeIntervalThres
            NSD(SI(I+1),chi)=false; % set coresponding to 0
        end
    end    
    
    %%% Apply accurate localization of peak  - by spline modeling
    if bAccuPeak
        SI=find(SD(:,chi));
        for k=1:length(SI)
            si=SI(k);
            binX=X(si-apSearchIntH:si+apSearchIntH,chi);
            bl=getbl3(binX);
            [~,idx]=max(bl);
            
            % assign new location
            SD(si,chi)=false;
            SD(si-apSearchIntH-1+idx, chi)=true;
        end
        
        if bNegPeak
        SI=find(NSD(:,chi));
        for k=1:length(SI)
            si=SI(k);
            binX=X(si-apSearchIntH:si+apSearchIntH,chi);
            bl=getbl3(binX);
            [~,idx]=min(bl);
            
            % assign new location
            NSD(si,chi)=false;
            NSD(si-apSearchIntH-1+idx, chi)=true;
        end
        end
    end
end


%%%%%%%%%%%%%% Spike windowed
sAmt=sum(SD);
A=cell(chAmt,1);
ptsAmt=preww+postww+1;
for chi=1:chAmt
    SI=find(SD(:,chi));
    if sAmt(chi)>0
        %%% Discard the spike too close to the edge of recording that window does not fit
        % * only need to consider the first and last one
        if SI(1)-preww<1
            SD(SI(1),chi)=false; SI(1)=[]; sAmt(chi)=sAmt(chi)-1;
        end
        if SI(end)+postww>pntAmt
            SD(SI(end),chi)=false; SI(end)=[]; sAmt(chi)=sAmt(chi)-1;
        end
        
        %%%
        % Output init
        A{chi}=zeros(ptsAmt,sAmt(chi));        
        for si=1:sAmt(chi)
            %%% data segmentation (different methods)
            dslow=SI(si)-preww; dshigh=SI(si)+postww;
            
%             % 1/N Y-align at maximum
%             temp=X(dslow:dshigh,chi);
%             temp=temp-X(I(si),chi);

            % 2/N Y-align at mean
            temp=X(dslow:dshigh,chi);
            temp=temp-mean(temp);
            
%             % 3/N XY-align at mean-crossing point
%             temp=temp-mean(temp);
%             temp=(temp<0);
%             idx=find(temp(I(si)-dslow:end),1);
%             dslow=I(si)+idx-fix(ww/2);
%             dshigh=I(si)+idx+fix(ww/2);
%             temp=X(dslow:dshigh,chi);
            
%             % 4/N XY-align
%             vmax=X(I(si),chi);
%             [vmin,idx]=min(X(I(si):I(si)+120,chi));
%             idx=I(si)+round(idx/2);
%             temp=X(idx-ww:idx+ww,chi);
%             temp=temp-(vmax+vmin)/2;
            
            %%% save data in output            
            A{chi}(:,si)=temp;
        end        
    end
end


%%%%%%%%%%%%%% Output the detection information
% 1. check the width of interval (longer than certain number can't be a real spike) *** not done yet; 
% 2. check the interval between spikes (but this should be considered
% with spike sorting, because 2 neighboring ones may belong to different
% cells !!!)
if bDQInfo
    DQ=cell(chAmt,1);
    spkInt=cell(chAmt,1);
    for chi=1:chAmt
        I=find(SD(:,chi));        
        
        %%% Width check        
        tQ=false(sAmt(chi),8);
        spkInt{chi}=zeros(sAmt,1); % store the intervals data, not complete
        for si=1:sAmt(chi)
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
                tQ(si,1)=true;
            end
        end
        
        %%% Spike interval check - assigned to the following spike which is
        %%% too close to previous one.
        for si=2:sAmt(chi)
            spkInt{chi}(si)=I(si)-I(si-1);
            if spkInt{chi}(si)<spikeIntervalThres
                tQ(si,2)=true;
            end
        end        
        
        DQ{chi}=tQ; % This works like: each row as a number
    end
    
    if nvarargout==1
        varargout{1}=DQ;
    end
end


%%%%%%%%%%%%%% Smoothing other peaks within the window
% * Method: replace original signal with its sqrt root - the bigger it was,
% the greater the reduction. At the same time ,remain the basic shape
if bSmoothOthers
    %%% Prepare length parameters to be used
    % post-peak window
    tp=round(postww*3/5);
    postw=[ptsAmt-tp+1,ptsAmt];
    % pre-peak window
    tp=round(preww/2);
    prew=[1,tp];
    % transition window half length
    thw=0.25/1000*srate;
    
    %%%
    for chi=1:chAmt        
        for si=2:sAmt(chi)
            if DQ{chi}(si,2) % if there is interval too small warning
                %%% If the interval (pts) is small enough. (1ms addtional
                % intended to add a residue to it.)
                if spkInt{chi}(si)<postww+0.001*srate
                    % First smooth the later half of previous spike.                  
                    temp=A{chi}(postw(1):postw(2),si-1);
                    I=(temp<0);
                    temp=sqrt(abs(temp));
                    temp(I)=-temp(I);
                    A{chi}(postw(1):postw(2),si-1)=temp;
                    % Transition smoothing                    
                    A{chi}(postw(1)-thw:postw(1)+thw,si-1)=smoo(A{chi}(postw(1)-thw:postw(1)+thw,si-1),2);
                    
                    % Then the first half of current spike
                    temp=A{chi}(prew(1):prew(2),si);
                    I=(temp<0);
                    temp=sqrt(abs(temp));
                    temp(I)=-temp(I);
                    A{chi}(prew(1):prew(2),si)=temp;
                    % Transition smoothing
                    A{chi}(prew(2)-thw:prew(2)+thw,si)=smoo(A{chi}(prew(2)-thw:prew(2)+thw,si),2);
                end                
            end
        end
    end
end
