% spike detection algorithm
%   [SD,A]=spikedetect(X,srate) 
% output the spike data SD (logical format{0,1} symbol), and windowed spike
% signals (for the later feature extraction)
%   [SD,A,SQ]=spikedetect(...), Q as quality marks
% SQ= detection quality information {channel}logical[spikeAmt,2]
% FALSE: normal; bit1 TRUE: too close to the previous spike; bit2 TRUE: over-threshold "plateau" too wide; 
% 4: low amplitude?
function [SD,SQ,varargout]=spike_detect(X,srate)
%%% Option parameter default setting - ordered according to procssing of data
% Whether to remove the outliner
useOLrm=true;
OLrmThres=20;% threshold for removal

% Whether to detect negative peaks.
bNegPeak=false;

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
    movLen=1; % length (s) of window   <<<<<<<<<<<<%%%%%%%%%%%%%
    movThresPos=4.2; % folds of estimated STD % 5 also works I think
    movThresNeg=5; % 
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
refraInterval=2.5; % (ms) 

% 是否用样条方法精确定位尖峰所在位置
bAccuPeak=true;
apSearchIntH=0.25; % half of the interval to search (left and right)

% window for spike waveform
preww=2.5; % pre-event window width (ms)
postww=5; % post-event

% 是否平滑主峰两侧的其他spike造成的峰
bSmoothOthers=true;

% Quality check related
spikeWidthThres=[0.1,3]; % (ms) should be in between
% 是否移除范围之外的spike
bBadWidthRemove=true;

spikeIntervalThres=6; % (ms) should be bigger than
% * spikeIntervalThres 应该比postww or preww大，这样才能在smoothing other spike 中提供保护
% （后者以DQ中spikeIntervalThres项目有记录为激活条件）


%%% Basic information
if nargout==3
    bAlign=true;
else
    bAlign=false;
end

[pntAmt,chAmt]=size(X);
movLen=srate*movLen;

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
[segm,sega]=cutseg([1,pntAmt],movLen);

movAvgSave=zeros(sega,chAmt);
for segi=1:sega
    % The RMS method need to separate the interval into smaller bins (20ms)
    if thresMethod==3        
        [sbseg,sba]=cutseg([1,segm(segi,2)-segm(segi,1)+1],srate*0.02);
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
        
        % Get the over-threshold part after processed average - all that matters
        % is the "movAvg" variable
        thresRes=(binX>movThresPos*ssd);        
        SD(segm(segi,1):segm(segi,2),chi)=thresRes;
        if bNegPeak
            thresRes=(binX<-movThresNeg*ssd);        
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

% Initiate The Spike data quality information， one cell for each channel.
% * It include 3 columns: 
% 1. marke polarity of spike: positive(1)/negative(0); 
% 2. check peak width (larger than certain number can't be a real spike): too wide(1)/normal(0);  
% 3. check the interval between spikes: too close to previous(1)/normal(0)
% (but this later should be considered with spike sorting, because 2 neighboring ones may belong to different cells)
SQ=cell(chAmt,1);

for chi=1:chAmt
    SQ{chi}=false(0,3);
    %%% Run through points
    spkcount=1;
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
            
            % Find the maximum signal location
            temp=X(ps:pe,chi);
            [~,idx]=max(temp);
            % set all points in this SD interval 0, leaving only the max loc
            SD(ps:pe,chi)=false;
            SD(ps+idx-1,chi)=true;
            
            % Polarity information in SQ (1.)
            SQ{chi}(spkcount,:)=false(1,3);
            SQ{chi}(spkcount,1)=true;
            
            % Peak width check
            if pe-ps>spikeWidthThres(2) || pe-ps<spikeWidthThres(1)
                SQ{chi}(spkcount,2)=true;
            end
            
            spkcount=spkcount+1;
        end
        % Negative
        if bNegPeak && NSD(pti,chi)            
            % mark the start edge loc
            ps=pti;
            % search for the end edge
            while pti<=pntAmt && NSD(pti,chi)
                pti=pti+1;
            end
            % mark the end edge
            pe=pti-1;

            % Find the maximum signal location
            temp=X(ps:pe,chi);
            [~,idx]=min(temp);
            % set all points in this NSD interval 0, leaving only the max loc
            NSD(ps:pe,chi)=false;
            NSD(ps+idx-1,chi)=true;
            
            % Polarity information in SQ (1.)
            SQ{chi}(spkcount,:)=false(1,3);
            
            % Peak width check
            if pe-ps>spikeWidthThres(2) || pe-ps<spikeWidthThres(1)
                SQ{chi}(spkcount,2)=true;
            end
            
            spkcount=spkcount+1;
        end
        pti=pti+1;
    end
    
    %%% Combine the Positive and Negative detections    
    if bNegPeak
        SD(:,chi)=SD(:,chi) | NSD(:,chi);
    end    
end
clear NSD


%%%%%%%%%%%%%%%%%%%% Checking - remove potentially bad spikes
rs=checkNegaInt*3; % range to calculate baseline for comparison
for chi=1:chAmt
    
    %%% Apply negative following check 
    if bCheckNegaFollow
        SI=find(SD(:,chi));
        rmlist=zeros(0,1);
        for k=1:length(SI)
            si=SI(k);
            % determine the segment of this point - in order to get the seg
            % baseline stored ahead.
%             Ilow=(si>=segm(:,1));
%             Ihigh=(si<=segm(:,2));
%             I=Ilow & Ihigh;
%             idx=find(I);            
%             bl=movAvgSave(idx,chi);
            
            temp=SD(si-rs:si+2*rs,chi);
            bl=mean(temp);
            
            % min value in the following period.
            [~,tp]=min(X(si:si+checkNegaInt,chi));
            if bl-X(si+tp,chi) > 2*(X(si,chi)-bl)
                SD(si,chi)=false;
                rmlist=[rmlist; k];
            end
        end        
        % negative peak here **************
        
        % Remove
        SQ{chi}(rmlist,:)=[];
    end
    
    %%% Apply statistical check
    if bCheckStatistic
        %
    end
    
    %%% 是否移除宽度在范围之外的spike
    if bBadWidthRemove
        SI=find(SD(:,chi));
        
        rmlist=find(SQ{chi}(:,2));
        
        SD(SI(rmlist),chi)=false;
        SQ{chi}(rmlist,:)=[];
    end

    %%% Apply refractory period
    % * not the real refractory, just too close so will affect
    % identification
    if useRefra
        SI=find(SD(:,chi));
        
        SIdif=diff(SI); % SIdif = interval
        rmlist=find(SIdif<refraInterval)+1; % find interval below spikeIntervalThres
        
        SD(SI(rmlist),chi)=false;
        SQ{chi}(rmlist,:)=[];
    end
        
    %%% Apply accurate localization of peak  - by spline modeling
    % * 这里需要考虑一下修正的波峰定位到另一个峰上的情况（比如两个峰很接近时），进而在“assign new location”时导致
    % 总spike检测数-1。要么通过强制性加入refractory，且refractory区间大于插值bl区间；要么采取特殊应对手段。
    if bAccuPeak
        SI=find(SD(:,chi));
        for k=1:length(SI)
            si=SI(k);
            binX=X(si-apSearchIntH:si+apSearchIntH,chi);
            bl=getbl3(binX);
            if SQ{chi}(k,1)
                [~,idx]=max(bl);
            else
                [~,idx]=min(bl);
            end
            
            % assign new location
            SD(si,chi)=false;
            SD(si-apSearchIntH-1+idx, chi)=true;
        end        
    end
end
sAmt=sum(SD);


%%%%%%%%%%%%%% Spike wave aligned to obtain window of data
if bAlign
A=cell(chAmt,1);
ptsAmt=preww+postww+1;
for chi=1:chAmt
    SI=find(SD(:,chi));
    if sAmt(chi)>0
        %%% Discard the spike too close to the edge of recording that window does not fit
        % * only need to consider the first and last one
        if SI(1)-preww<1
            SD(SI(1),chi)=false; SI(1)=[]; SQ{chi}(1,:)=[]; sAmt(chi)=sAmt(chi)-1;
        end
        if SI(end)+postww>pntAmt
            SD(SI(end),chi)=false; SI(end)=[]; SQ{chi}(end,:)=[]; sAmt(chi)=sAmt(chi)-1;
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
%     %%%%### for debugging
%     if sum(SD(:,chi))~=size(SQ{chi},1)
%         error('no match');
%     end
end
end

%%%%%%%%%%%%%% Smoothing other peaks within the window
% * Method: replace original signal with its sqrt root - the bigger it was,
% the greater the reduction. At the same time ,remain the basic shape

%%% First add the ISI check to the SQ - assigned to the following spike which is
% too close to previous one - disregard of positive/negative.
spkInt=cell(chAmt,1);
for chi=1:chAmt
    % Get interval
    I=find(SD(:,chi));
    spkInt{chi}=[0;diff(I)];
    % Spike interval check
    I=(spkInt{chi}<spikeIntervalThres);
    SQ{chi}(I,3)=true;
%     for si=2:sAmt(chi)
%         spkInt{chi}(si)=I(si)-I(si-1);
%         if spkInt{chi}(si)<spikeIntervalThres
%             SQ{chi}(si,3)=true;
%         end
%     end
end

%%%
if bSmoothOthers && bAlign
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
            if SQ{chi}(si,2) % if there is interval too small warning
                %%% If the interval (pts) is small enough. (addtional 1ms 
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

if bAlign, varargout{1}=A; end
