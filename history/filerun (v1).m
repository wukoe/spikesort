% "run-through" for operating the file piece by piece.
%   filerun(procop,outfn,infn,paras): method, output, input, parameters
%  input: original file; output: detected spike trains
% MAJORLY TRY A CHANNEL BASED SIZE STRATEGY
function success=filerun(procop,outfn,varargin)
%%%%%%%% Data Initiation
success=false;
% Default parameter setting
bOverWrite=true;

% Load the setting information
% find the location of spikesort function
sspath=which('spikesort');
if isempty(sspath),    error('can not locate the file for parameters'); end
sspath=fname(sspath); 
sspath=[sspath,'rsdata.mat'];
% load parameters
load(sspath); % usually contains "paras", "info"
bDirectMem=paras.bServer; % whether direct memory instead of using "matfile"

% locate the input file one by one
nvarargin=length(varargin);
for k=1:nvarargin
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

% locate the output file - in and out file use the same location
nstr=[outfn,'.mat'];
if exist(nstr,'file')
    if bOverWrite,   delete(nstr); end
end
if bDirectMem
    outfile=struct();
else
    outfile=matfile(nstr,'Writable',true);
end


%%%%%%%%%%% The preprocessing part - !!! The separation of initiation and
% main processing is intended to fit well with both channel-based and
% segment-based processing (although the later is sealed here).
% * this part also check the input option.
switch procop
    case 'filt'
        ptsAmt=info.ptsAmt;chAmt=info.chAmt;
        outfile.srate=infile1.srate;
        outfile.T=infile1.T;
        outfile.X=zeros(ptsAmt,1);
        
    case 'spike detect'
%         outfile.T=infile1.T;
        ptsAmt=info.ptsAmt;chAmt=info.chAmt;
        outfile.SD=false(ptsAmt,1);
        outfile.SQ=cell(chAmt,1);        
        
    case 'spike align'
        chAmt=info.chAmt;
        outfile.A=cell(chAmt,1);
        %>>> infile1=X; infile2=SD;
        
    case 'spike feature'
        chAmt=info.chAmt;
        outfile.SF=cell(chAmt,1);
        
    case 'cluster'
        ptsAmt=info.ptsAmt; chAmt=info.chAmt;
        T=infile2.T;
        if length(T)~=ptsAmt, error('may be loading wrong data'); end
        
        totalTime=(T(ptsAmt,1)-T(1,1))/1000; % total time span of data        
        
        outfile.time=[T(1,1),T(ptsAmt,1)];
        outfile.CSI=cell(chAmt,1); % Channel SI: record which of cluster identity of each spike in original channel.
        % note there may be 0 in CSI - those not assigned to any cluster.
        
%     case 'final'
%         tp=size(infile1,'SSD');    neuAmt=tp(2);
%         T=infile2.T;
%         
%         outfile.NSD=cell(neuAmt,1);
%         outfile.neuCh=infile1.neuCh;
%         outfile.time=[T(1),T(end)];        
        
    otherwise
        error('invalid option');
end


%%%%%%%%%% Main Processing
switch procop    
    case 'filt'
        disp('filt >>>');
        chID=zeros(chAmt,1); % the channel ID 
        outchi=1;        
        for chi=1:chAmt
            temp=filtfilt(paras.fb,paras.fa,double(infile1.X(:,chi)));
            
%             removing low-freq noise & base line of signal %%%%%%% Use
%             that method, and filtfilt will focus on high freq later
            
            % When the noise level (measured by STD) overpass a threshold
            if paras.bRmc && std(temp)>paras.rmChThres 
                fprintf('X');
            else
                outfile.X(:,outchi)=temp;
                chID(outchi)=info.chID(chi);
                outchi=outchi+1;
                fprintf('|');
            end
        end        
        fprintf('\n');
        % remove extra space of chID if any
        if outchi<chAmt,    chID(outchi+1:chAmt)=[];        end
        outfile.chID=chID;
        
    case 'spike detect'
        disp('spike detect >>>');         
        for chi=1:chAmt
            [outfile.SD(:,chi),temp,outfile.SQ(chi,1)]=spike_detect(infile1.X(:,chi),info.srate); 
            outfile.A(chi,1)=temp;
            fprintf('|');
        end
        fprintf('\n');
        outfile.T=infile1.T;
        outfile.chID=infile1.chID;
        
    case 'spike align'
        disp('spike align >>>');
        for chi=1:chAmt
            outfile.A(chi,1)=spike_align(infile1.X(:,chi),infile2.SD(:,chi),info.srate);
            fprintf('|');
        end
        fprintf('\n');
        
    case 'spike feature'
        disp('spike feature >>>');
        for chi=1:chAmt
            temp=infile1.A(chi,1);
            temp{1}=spike_feature(temp{1});
            outfile.SF(chi,1)=temp;
            fprintf('|');
        end
        fprintf('\n');
        
    case 'cluster'
        disp('spike clustering >>>');
        % neuCh (neuron channel) is the table noting from which of the
        % original electrodes the neuron is recorded.
        neuCh=zeros(1,1);
        % NSD is neuron spk data (index format)
        NSD=cell(1,1);
        SI=cell(1,1);
        
        count=1;
        for chi=1:chAmt
            temp=infile1.SF(chi,1);
            spkTime=T(infile2.SD(:,chi));
            spkData=infile3.A(chi,1);
            [SI{1},lbinfo]=spike_cluster(temp{1},'kmeans',totalTime,spkTime,spkData{1});
             % it is possible that non is left; the "neuCh" is based on
             % counting.
            if lbinfo.cAmt>0
                spLoc=T(infile2.SD(:,chi));            
                for k=1:lbinfo.cAmt
                    % 1/2 logical format
%                     I=spLoc(lbinfo.ids{k}); % current neuron's activity index of the electrode spike event.
%                     temp=false(ptsAmt,1);  temp(I)=true;
%                     outfile.SSD(:,count)=temp;
                    % 2/2 index format
                    % sort the index in ascending order, since lbinfo.ids{} does not guarantee monotonics.
                    temp=sort(lbinfo.ids{k}); 
                    NSD{count,1}=spLoc(temp);
                    neuCh(count,1)=infile2.chID(chi,1);

                    count=count+1;
                end
            end
            
            outfile.CSI(chi,1)=SI;
            
            fprintf('|');
        end
        fprintf('\n');
        
        outfile.NSD=NSD;
        outfile.neuCh=neuCh;        
        
%     case 'final'
%         disp('final, almost done >>>');
%         
%         % transform into the spike time index format
%         for ni=1:neuAmt
%             temp{1}=T(infile1.SSD(:,ni));
%             outfile.NSD(ni,1)=temp;
%         end
end

%%%%%%%%%%%%% Save the output data if runs in direct memory mode
if bDirectMem
    save(nstr,'-v7.3','-struct','outfile');
end
disp('>>> results saved, current running finished');

% % %%%%%%%%%%% 此部分用于分段处理，封存！
% % % to choose specific start and end point
% % startLine=floor(startLoc/100 * lineAmt);
% % if startLine+readLineNum>=lineAmt
% %     readLineNum=lineAmt-startLine;
% % end
% % tp=startLine * lineLen + fhe;
% % fseek(fid,tp,-1);
% % fgetl(fid);

% % % calculate how many segment is needed
% % [bSeg,bAmt]=cutseg(ptsAmt,blockSize);
% % fprintf('total number of pieces: %d\n',bAmt);
% % for bi=1:bAmt
% %     ptilow=bSeg(bi,1); ptihigh=bSeg(bi,2);
% %     
% %     
% %     
% % end
% % 
% % %%%%%%%%% Evaluate the block size
% % % Get system spare memory quantity
% % try
% %     [~,mem]=memory();
% %     mem=mem.PhysicalMemory.Available;
% % catch
% %     % If system mem info is unavailble, suppose it is 1GB on PC, 5GB on
% %     % server
% %     disp('warning: system memory information unavailable.');
% %     if bServer
% %         mem=10*1024^3; % won't be much when multiplied by the ratio
% %     else
% %         mem=1*1024^3;
% %     end
% % end
% % fprintf(fid,'memory to be used: %0.2g GB\n',mem*blockSize/1024^3);
% % 
% % blockSize=floor(mem*blockSize/info.lineLen);
% % fprintf(fid,'The processing block size: %d\n',blockSize);
