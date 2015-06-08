% Delete the duplicates from conduction signal
%   [newST,marklb,rmlb]=spikemerge(ST,SA,info)
% options:
% 'ext'=0.0005 (s)
% 'auto thres'=true
function [ST,marklb,rmlb]=spikemerge(ST,SA,info,varargin)
% Standard to find cluster sequence.
ext=0.0005; %0.5/0.25ms
% optSelect=2; % way to determine the representative of the cluster
bAutoThres=true; % whether use algorithm-determined adaptive threshold for pair number.
minThres=4/60*info.TimeSpan; % static threshold if chosen to use.
flagUserTimeSpan=false;
bRmInt=false;

if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'ext'
                ext=pinfo{parai};
            case 'time span'
                timeSpan=pinfo{parai};
                flagUserTimeSpan=true;
            case 'auto thres'
                bAutoThres=pinfo{parai};
            otherwise
                error('unidentified options');
        end
    end
end

%%% Process
chAmt=length(ST);
% Get number of firing in each channel
sAmt=cellstat(ST,'length');
% Check number of SD SA
assert(isequal(sAmt,cellstat(SA,'length')),'merging: SD and SA number not equal');

if ~flagUserTimeSpan
    if bRmInt        
        % First remove the large "empty" segments of firing segments
        % (most typical case: the inter-burst interval)
        for chi=1:chAmt
            bl=burstlet(ST{chi});
        end
    else
        timeSpan=max(cellstat(ST,'max'))-min(cellstat(ST,'min'));
    end
end

%%% Identify the spikes involved in conduction signal.
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
    if ~bAutoThres
        if length(ST{chi})<minThres
            continue
        end
    else
        if length(ST{chi})<3
            continue
        end
    end
    
    % Find all tight-following spikes (under ext) in the channel pairs
    [Ichi,Isd,rnum]=cspair(ST{chi},ST(chi+1:end),ext); % ����Ҫ��ȫ�����ɨ�裬ֻ��Ҫ�Ժ���ĵ�Ԫ���м��ɡ���Ϊ��ǰ����Ѿ����й��ˡ�
    
    if bAutoThres
        pthr=1e-6;
        
        sdnum=length(Isd);
        ts=sAmt(chi+1:end);
        foI=false(sdnum,1);
        for m=1:sdnum            
            [P,expectnum]=probable_rs(timeSpan,[sAmt(chi),ts(m)],ext,rnum(m));
            % * when 2ch have spike 10000&1000, expect~30; when have
            % 10000&10000, expect~300-400.
            if isnan(P)
                error('check this out: P calculation numeric problem');
            end
            if rnum(m)>expectnum && P<pthr 
                % ��һ��������Ϊ�˱�֤P�㹻С��ԭ������Ϊrnumλ�ڸ��ʷֲ���ֵ����һ�ࣨ������󳬹����ģ���µ�Ԥ�ڣ���
                % �����෴������ԶС��Ԥ�ڣ���
                foI(m)=true;
            end            
        end
        foI=find(foI); foAmt=length(foI);
    else
        foI=find(rnum>=minThres); foAmt=sum(rnum>=minThres);        
    end
    if foAmt==0 % ��û���κδﵽ�㹻������conductionͨ�����������channel.
        continue
    end
    
    % Leave only effective channels.
    Ichi=Ichi(:,foI); Isd=Isd(foI); rnum=rnum(foI);
    % ����ÿ��ͨ���ġ�following��spike���Ҹ��Եķ�ֵ��
    ampchi=zeros(foAmt,1);
    ampsd=ampchi;
    for fi=1:foAmt
        % itself
        ampchi(fi)=mean(abs(SA{chi}(Ichi(:,fi))));
        % conducters
        ampsd(fi)=mean(abs(SA{foI(fi)+chi}(Isd{fi})));
    end    
%     % �������ampchi���ı��Ӧ��ampsd��foI
%     [ampchi,ampI]=sort(ampchi,'descend');
%     ampsd=ampsd(ampI);
%     foI=foI(ampI); 
%     Ichi=Ichi(:,ampI); Isd=Isd(ampI); rnum=rnum(ampI);
    for fi=1:foAmt
        % Comparing rule
        % if the amplitude difference is larger than 1.25 times of the
        % other, select by larger amplitude; else if the amplitude
        % difference is not large enough, select by more total spike number
        % (also only suitable for sorted data).
        % * ����ȫ����chi�б�������spikeȺ�����ѷ���/δ����ı��������Ƿ���Щ�ҵ�����Ϊ������Ⱥ�壬
        % ����detectionȱ�ݵ���ս��
        if ampchi(fi)>ampsd(fi)*1.25 || ampsd(fi)>ampchi(fi)*1.25
            if ampchi(fi)>ampsd(fi)
                flagKeepchi=true;
            else
                flagKeepchi=false;
            end
        else
            %
            if sAmt(chi)>sAmt(foI(fi)+chi)
                flagKeepchi=true;
            else
                flagKeepchi=false;
            end
        end
        
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

%%% Delete from the spike train & label of representative identity marker
for chi=1:chAmt
    ST{chi}(rmlb{chi})=[];
    marklb{chi}(rmlb{chi})=0;
end

end % of main

%     if optSelect==1 % keep channel with maximum total firing.
%         % percentage of repeats to each unit's total firing number
%         P=stat.seqcount(seqi)./sAmt(seq);
%         [~,idx]=max(P);
%     else % keep channel with maximum spike amplitude.
%         % Get amplitude of spikes in each channel (use only spikes participating current cluster,
%         % remove possible artifact)
%         sAmp=zeros(seqlen,1);
%         for k=1:seqlen
%             I=seqCS{seqi}(:,k);
%             temp=abs(SA{seq(k)}(I));
%             I=outlier_detect(temp);
%             temp=temp(~I);
%             sAmp(k)=mean(temp);
%         end
%         [~,idx]=max(sAmp);
%     end