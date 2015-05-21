% Delete the duplicates from conduction signal
%   [newST,marklb,rmlb]=spikemerge(ST,info,SA)
% options:
% 'ext'=0.0005 (s)
% 'auto thres'=true
function [ST,marklb,rmlb]=spikemerge(ST,SA,info,varargin)
% Standard to find cluster sequence.
ext=0.0005; %0.5/0.25ms
% optSelect=2; % way to determine the representative of the cluster
bAutoThres=false; % whether use algorithm-determined adaptive threshold for pair number.
minThres=4/60*info.TimeSpan; % static threshold if chosen to use.

if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'ext'
                ext=pinfo{parai};
%             case 'total'
%                 totaltimelimit=pinfo{parai};
            case 'auto thres'
                bAutoThres=pinfo{parai};
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
    end
    
    % Find all tight-following spikes (under ext) in the channel pairs
    [Ichi,Isd,rnum]=cspair(ST{chi},ST(chi+1:end),ext); % ����Ҫ��ȫ�����ɨ�裬ֻ��Ҫ�Ժ���ĵ�Ԫ���м��ɡ���Ϊ��ǰ����Ѿ����й��ˡ�
    
    if bAutoThres
        pthr=1e-6;
        
        sdnum=length(Isd);
        ts=sAmt(chi+1:end);
        foI=false(sdnum,1);
        for m=1:sdnum
            P=probable_rs(info.TimeSpan,[sAmt(chi),ts(m)],ext,rnum(m));
            % * when 2ch have spike 10000&1000, expect~30; when have
            % 10000&10000, expect~300-400.
            if P<pthr
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
    
    % ����ÿ��ͨ���ġ�following��spike���Ҹ��Եķ�ֵ��
    ampchi=zeros(foAmt,1);
    ampsd=ampchi;
    for fi=1:foAmt
        idx=foI(fi);
        % itself
        ampchi(fi)=mean(abs(SA{chi}(Ichi(:,idx))));
        % conducters
        ampsd(fi)=mean(abs(SA{idx+chi}(Isd{idx})));
    end    
    % �������ampchi���ı��Ӧ��ampsd��foI
    [ampchi,ampI]=sort(ampchi,'descend');
    ampsd=ampsd(ampI);
    foI=foI(ampI); 
    Ichi=Ichi(:,foI); Isd=Isd(foI); rnum=rnum(foI);    
    for fi=1:foAmt
        if ampchi(fi)>ampsd(fi) %����chi            
            % ����������ֵ�Ŀ���Ǹ���chi�б�������spikeȺ�����ѷ���/δ����ı��������Ƿ���Щ�ҵ�����Ϊ������Ⱥ�塣
            % ����һ�뷨����detectionȱ�ݵ���ս��Ŀǰ��ʱ���ӡ�
%             tp=cslb{chi}(Ichi(:,fi)); %������conduction��cslb
%             nU=sum(tp==0); % δ�����spike����
%             % ��������ȡ�����ѷ�����������������cslb��spike������
%             tp(tp==0)=[];
%             lb=reabylb(tp);
%             [nE,idx]=max(lb.typeAmt);
%             existlb=lb.types(idx);
%             % ��ȡexistlb��������chiͨ��spike���������
%             nEall=sum(cslb{chi}==existlb);            
%             % ��������������������new lb
%             if nE>nEall*0.9 % ��ȫ��nE��Ϊ
%             end

            % mark the "marker" spikes
            marklb{chi}(Ichi(:,fi))=1;            
            % ���Ҫɾ����sd(fi)
            rmlb{foI(fi)+chi}(Isd{fi})=true;
        else
            rmlb{chi}(Ichi(:,fi))=true;            
            marklb{foI(fi)+chi}(Isd{fi})=1;
        end
    end
    fprintf('|');
end
fprintf('\n');

%%% Delete from the spike train & representative identity mark lb
for chi=1:chAmt
    ST{chi}(rmlb{chi})=[];
    marklb{chi}(rmlb{chi})=0;
end

end % of main