% Delete the duplicates from conduction signal
%   [newST,marklb,rmlb]=spikemerge(ST,SA,info)
% options:
% 'ext'=0.0005 (s)
% 'auto thres'=true
function [ST,marklb,rmlb]=spikemerge(ST,SA,info,varargin)
% Standard to find cluster sequence.
ext=0.0008; %0.5/0.25ms
extVarThre=0.0001; 
% optSelect=2; % way to determine the representative of the cluster
param=struct();
param.bAutoThres=true; % whether use algorithm-determined adaptive threshold for pair number.
param.numPthr=1e-6; % for with auto threshold
minThres=4/60*info.TimeSpan; % static threshold if chosen to use.
flagUserTimeSpan=false;
bRmIntv=false; % whether remove the large "empty" segments of spike train
% spike time difference jitter limit.
param.DTJthr=0.15/1000;%(ms/1000) = 3 points at 20000HZ SR. <<<

% User
if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'ext'
                ext=pinfo{parai};
            case 'time span'
                param.timeSpan=pinfo{parai};
                flagUserTimeSpan=true;
            case 'auto thres'
                param.bAutoThres=pinfo{parai};
            otherwise
                error('unidentified options');
        end
    end
end

% Process
chAmt=length(ST);
% Get number of firing in each channel
sAmt=cellstat(ST,'length');
% Check number of SD SA
assert(isequal(sAmt,cellstat(SA,'length')),'merging: SD and SA number not equal');

if ~flagUserTimeSpan
    if bRmIntv
        % First remove the large "empty" segments of spike train
        % (most typical case: the inter-burst interval)
        for chi=1:chAmt
            % <<<
        end
    else
        param.timeSpan=max(cellstat(ST,'max'))-min(cellstat(ST,'min'));
    end
end


%%%%%%%%%%%%% Identify the spikes involved in conduction signal.
% Initialize the markers (Spike remain/delete, etc).
% cslb=cell(chAmt,1);% �������ͨ����
marklb=cell(chAmt,1);% ����spike�Ƿ�ѡΪ����ı��
rmlb=cell(chAmt,1); % ����spike�Ƿ�ɾ��
for chi=1:chAmt
%     cslb{chi}=zeros(sAmt(chi),1);
    marklb{chi}=zeros(sAmt(chi),1);
    rmlb{chi}=false(sAmt(chi),1);
end
for chi=1:chAmt
    if ~param.bAutoThres
        if length(ST{chi})<minThres
            continue
        end
    else
        if length(ST{chi})<3
            continue
        end
    end
    
    %%% Find all tight-following spikes (under ext) in the channel pairs
    [Ichi,Isd,~,foI]=cspair(ST{chi},ST(chi+1:end),ext,param); % ����Ҫ��ȫ�����ɨ�裬ֻ��Ҫ�Ժ���ĵ�Ԫ���м��ɡ���Ϊ��ǰ����Ѿ����й��ˡ�
    
%     %%% Number threshold.
%     if bAutoThres
%         sdnum=length(Isd);
% %         ts=sAmt(chi+1:end);
%         foI=false(sdnum,1);% indicator of "follower"
%         for m=1:sdnum
%             [P,expectnum]=probable_rs(timeSpan,[sAmt(chi),sAmt(chi+m)],ext,rnum(m));
%             % * when 2ch have spike 10000&1000, expect~30; when have
%             % 10000&10000, expect~300-400.
%             if rnum(m)>expectnum && P<numPthr
%                 % ��һ��������Ϊ�˱�֤����һ����ʵ��P�㹻С��ԭ������Ϊrnumλ�ڸ��ʷֲ���ֵ���Ҳࣨ������󳬹����ģ���µ�Ԥ�ڣ���
%                 % �����෴������ԶС��Ԥ�ڣ���
%                 foI(m)=true;
%             end
%         end
%         foI=find(foI); foAmt=length(foI);
%     else
%         foI=find(rnum>=minThres); foAmt=sum(rnum>=minThres);
%     end
%     % ��û���κ�����������conductionͨ�������channel(chi)��������ص�CS������.
%     if foAmt==0 
%         continue
%     end
%     % Keep only effective CS channels's information.
%     Ichi=Ichi(:,foI); Isd=Isd(foI); rnum=rnum(foI);
%     
%     %%% Time difference threshold
%     for fi=1:foAmt
%         % Time difference between each pair.
%         dt=ST{chi}(Ichi(:,fi))-ST{foI(fi)+chi}(Isd{fi});
%         % Get mean
%         dtm=mean(dt);
%         % Delete any pair whose difference deviate from this mean more than
%         % threshold.
%         I=(abs(dt-dtm)>DTJthr);
%         tp=find(Ichi(:,fi));
%         Ichi(tp(I),fi)=false;
%         % * removing Isd must be placed later ?
%         tp=find(Isd{fi});
%         Isd{fi}(tp(I))=false;
%         % Update repeat number.
%         rnum(fi)=sum(Ichi(:,fi));
%     end
% %     % Remove the unqualified spikes from Isd list
% %     for fi=1:foAmt
% %         tp=find(Isd{fi});
% %         Isd{fi}(tp(I))=false;
% %     end
%     
%     % Check for second time whether the remaining still fullfill number limit.
%     if bAutoThres
%         seleI=false(foAmt,1);
%         for fi=1:foAmt
%             [P,expectnum]=probable_rs(timeSpan,[sAmt(chi),sAmt(chi+foI(fi))],ext,rnum(fi));
%             if rnum(fi)>expectnum && P<numPthr
%                 seleI(fi)=true;
%             end
%         end
%     else
%         seleI=(rnum>=minThres);
%     end
%     foI=foI(seleI); foAmt=sum(seleI);
%     % Keep only effective CS channels's information.
%     Ichi=Ichi(:,seleI); Isd=Isd(seleI); rnum=rnum(seleI);
%     % �ٴΣ���û���κ�����������conductionͨ�����������channel(chi).
    foAmt=length(foI);
    if foAmt==0 
        continue
    end    
    
    %%% Determine which to delete of each pair.
    % * Comparing rule:
    % if the amplitude difference is larger than 1.25 times of the
    % other, select by larger amplitude; else if the amplitude
    % difference is not large enough, select by more total spike number
    % (also only suitable for sorted data).
    % * ����ȫ����chi�б�������spikeȺ�����ѷ���/δ����ı��������Ƿ���Щ�ҵ�����Ϊ������Ⱥ�壬
    % ����detectionȱ�ݵ���ս��
    
    % ����ÿ��ͨ���ġ�following��spike���Ҹ��Եķ�ֵ��
    ampchi=zeros(foAmt,1);
    ampsd=ampchi;
    for fi=1:foAmt
        % itself
        ampchi(fi)=mean(abs(SA{chi}(Ichi(:,fi))));
        % conducters
        ampsd(fi)=mean(abs(SA{foI(fi)+chi}(Isd{fi})));
    end
    
    % Apply the comparing rule.
    for fi=1:foAmt
%         if ampchi(fi)>ampsd(fi)*1.25 || ampsd(fi)>ampchi(fi)*1.25
            if ampchi(fi)>ampsd(fi)
                flagKeepchi=true;
            else
                flagKeepchi=false;
            end
%         else
%             %
%             if sAmt(chi)>sAmt(foI(fi)+chi)
%                 flagKeepchi=true;
%             else
%                 flagKeepchi=false;
%             end
%         end
        
        if  flagKeepchi %����chi
            % mark the "marker" spikes
            marklb{chi}(Ichi(:,fi))=1;
            % ���Ҫɾ����sd(fi)
            rmlb{foI(fi)+chi}(Isd{fi})=true;
        else
            marklb{foI(fi)+chi}(Isd{fi})=1;
            rmlb{chi}(Ichi(:,fi))=true;
        end
    end
    fprintf('|');
end
fprintf('\n');

% Delete from the spike train & label of representative identity marker
for chi=1:chAmt
    ST{chi}(rmlb{chi})=[];
    marklb{chi}(rmlb{chi})=0;
end

end % of main