% "run-through" for operating the file piece by piece.
%   filerun(procop,outfn,infn,paras): method, output, input, parameters
%  input: original file; output: detected spike trains
% MAJORLY TRY A CHANNEL BASED SIZE STRATEGY
function success=filerun(procop,outfn,varargin)
%%%%%%%%%% Initiation
%%% Default
success=false;
% Default parameter setting
bOverWrite=true;
bUpdateStat=false; % whether to update the 'runstat.mat'

%%% Load data information and parameters
% Locate the run state file
sspath=which([varargin{1},'.mat']); % same as spikesort function
if isempty(sspath),    error('can not locate the file for parameters'); end
sspath=fname(sspath);
sspath=[sspath,'rsdata.mat'];
% loading
load(sspath); % usually contains "paras", "info"
bDirectMem=paras.bDirectMem; % whether direct memory instead of using "matfile"

%%% Locate the input data (in file format) one by one
nvarargin=length(varargin);
for k=1:nvarargin
    if ischar(varargin{k});
    nstr=[varargin{k},'.mat'];    
    if exist(nstr,'file')
        if bDirectMem
            tstr=sprintf('infile%d=load(nstr);',k);% infile1, infile2, ... according to 
        else
            tstr=sprintf('infile%d=matfile(nstr);',k);
        end
        eval(tstr);
    else
        error('did not find #%d input file!',k);
    end
    end
end

%%% Locate the output file - in and out file use the same location
nstr=[outfn,'.mat'];
if exist(nstr,'file')
    if bOverWrite,   delete(nstr); end
end
if bDirectMem
    outfile=struct();
else
    outfile=matfile(nstr,'Writable',true);
end


%%%%%%%%%%% The initialization part - 
% !!! The separation of initiation and main processing 
% is intended to fit well with both channel-based and
% segment-based processing (although the later is sealed here).
% * this part also check the input option.
switch procop
    case 'filt'
        chAmt=info.chAmt;
        outfile.srate=infile1.srate;
        outfile.T=infile1.T;
        outfile.X=zeros(info.ptsAmt,chAmt);
        
    case 'spike detect'
        chAmt=info.chAmt;
        outfile.SD=cell(chAmt,1);
        outfile.SQ=cell(chAmt,1);        
        
    case 'spike feature'
        chAmt=info.chAmt;
        outfile.A=cell(chAmt,1);
        outfile.SF=cell(chAmt,1);
        
    case 'cluster'
        chAmt=info.chAmt;
        T=varargin{3};
        if length(T)~=info.ptsAmt, error('may be loading wrong data'); end
        
        totalTime=(T(info.ptsAmt,1)-T(1,1))/1000; % total time span of data        
        
        outfile.time=[T(1,1),T(info.ptsAmt,1)];
        outfile.CSI=cell(chAmt,1); % Channel SI: record which of cluster identity of each spike in original channel.
        % note there may be 0 in CSI - those not assigned to any cluster.
        
    otherwise
        error('invalid option');
end


%%%%%%%%%% Main Processing
switch procop
    %%%
    case 'filt'        
        rmMark=false(chAmt,1);
        for chi=1:chAmt
            temp=filtfilt(paras.fb,paras.fa,double(infile1.X(:,chi)));

            % When the noise level (measured by STD) overpass a threshold
            if paras.bRmc && std(temp)>paras.rmChThres 
                fprintf('X');
                rmMark(chi)=true;
            else
                outfile.X(:,chi)=temp;                
                fprintf('|');
            end
        end        
        fprintf('\n');
        % remove extra space of chLabel if any
        temp=info.chLabel; temp(rmMark)=[];
        outfile.chLabel=temp;% the channel ID 
        info.newChLabel=info.chLabel;
       
        
    %%%
    case 'spike detect'
        disp('spike detect >>>');
        temp=cell(2,1);
        for chi=1:chAmt
            [temp{1},temp{2}]=spike_detect(infile1.X(:,chi),info.srate);
            outfile.SD(chi,1)=temp(1); outfile.SQ(chi,1)=temp(2);
            fprintf('|');
        end
        fprintf('\n');
        
        %%% Separat the positive and negative peaks to separate channels
        %%% here. ¡¾need¡¿
        if paras.bPNSep
            disp('update channels >>>');
            newSD=cell(0,1);%false(ptsAmt,0);
            newA=cell(0,1);
            newSQ=cell(0,1);
            newchLabel=zeros(0,1);

            chcount=0;
            for chi=1:chAmt
                PN=outfile.SQ(chi,1);
                PN=PN{1}(:,1); % now [snum*1]: only the posi/nega mark
                % position of all the spikes
                SI=outfile.SD{chi,1};

                pnum=sum(PN); % number of positive peaks
                if pnum>0
                    chcount=chcount+1;

                    % assign to the newSD
                    newSD{chcount,1}=SI(PN);
                    % assign to new SQ
                    temp=outfile.SQ(chi,1);
                    newSQ{chcount,1}=temp{1}(PN,:);
                    % assign new chLabel
                    newchLabel(chcount,1)=infile1.chLabel(chi,1);                
                end

                if pnum<length(SI) % number of negative peaks
                    chcount=chcount+1;

                    PN=~PN;
                    % assign to the newSD
                    newSD{chcount,1}=SI(PN);
                    % assign to new SQ
                    temp=outfile.SQ(chi,1);
                    newSQ{chcount,1}=temp{1}(PN,:);
                    % assign new chLabel
                    newchLabel(chcount,1)=infile1.chLabel(chi,1);                
                end
            end

            outfile.SD=newSD;
            outfile.A=newA;
            outfile.SQ=newSQ;
            outfile.chLabel=newchLabel;

            info.chAmt=chcount;
            info.newChLabel=newchLabel;
            bUpdateStat=true;
        end
        
       
    %%% Get features
    case 'spike feature'
        disp('spike feature >>>');
        for chi=1:chAmt
            temp=infile1.A(chi,1);
            if isempty(temp{1})
                outfile.SF(chi,1)=temp;
            else
                temp{1}=spike_feature(temp{1});
                outfile.SF(chi,1)=temp;
            end
            fprintf('|');
        end
        fprintf('\n');
        
       
    %%%
    case 'cluster'
        disp('spike clustering >>>');
        % neuCh (neuron channel) is the table noting from which of the
        % original electrodes the neuron is recorded.
        neuCh=zeros(1,1);
        % NSD is neuron spk data (index format)
        NSD=cell(1,1);        
        
        count=1;
        for chi=1:chAmt
            temp=infile1.SF(chi,1);
            spkTime=T(infile2.SD{chi});
            spkData=infile1.A(chi,1);
            if length(spkTime)<=2
                NSD{count,1}=spkTime;
                outfile.CSI(chi,1)={ones(length(spkTime),1)};
                neuCh(count,1)=chi;
                
                count=count+1;
            else
            [SI,lbinfo]=spike_cluster(temp{1},'kmeans',totalTime,spkTime,spkData{1});
             % it is possible that non is left; the "neuCh" is based on
             % counting.
            if lbinfo.cAmt>0
%                 spLoc=find(infile2.SD(:,chi));            
                for k=1:lbinfo.cAmt
                    % 1/2 logical format
%                     I=spLoc(lbinfo.ids{k}); % current neuron's activity index of the electrode spike event.
%                     temp=false(ptsAmt,1);  temp(I)=true;
%                     outfile.SSD(:,count)=temp;
                    % 2/2 Time index format
                    % sort the index in ascending order, since lbinfo.ids{} does not guarantee monotonics.
                    temp=sort(lbinfo.ids{k}); 
                    NSD{count,1}=spkTime(temp); %T(spLoc(temp)); %%%%%%%%%?
                    neuCh(count,1)=chi; %info.chLabel{chi};

                    count=count+1;
                end
            end
            
            outfile.CSI(chi,1)={SI};
            end
            
            fprintf('|');
        end
        fprintf('\n');
        
        outfile.NSD=NSD;
        outfile.neuCh=neuCh;
        outfile.T=T;        
end

%%%%%%%%%%%%% Final closing
%%% Save the output data if runs in direct memory mode
if bDirectMem
    save(nstr,'-v7.3','-struct','outfile');
end

%%% Update the data information if necessary
if bUpdateStat
    save(sspath,'-append','info');
end
disp('>>> results saved, current running finished');
