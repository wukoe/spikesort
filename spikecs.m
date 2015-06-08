% Locate conduction signal in the original spike train data.
% Use results from spikemerge().
function [seq,seqMark,seqCS,stat]=spikecs(ST,marklb,info)
chAmt=length(ST);
ext=0.001;
bAutoThres=true;
minThres=4/60*info.TimeSpan;

%%% Find the trace by the marked "top" units
stat=struct();
% Get number of "marker" spikes of each unit.
markSDnum=zeros(chAmt,1);
for chi=1:chAmt
    markSDnum(chi)=sum(marklb{chi}==1);
end
% Filt out units with enough marker.
markI=find(markSDnum>0); markAmt=length(markI);

% End game if no marker neuron is found.
if markAmt==0
    seq=[];    seqMark=[];    seqCS=[];
    stat.marklb=marklb;
    return
end

%%% Use the markers to re-scan the SD.
nsd2=ST(markI);
marklb=marklb(markI); stat.markSDnum=markSDnum(markI);
seqMark=markI;% marker neuron ID.

sAmt=cellstat(ST,'length');

% Take mark neurons in turn to scan.
seq=cell(markAmt,1);
seqCS=cell(markAmt,1);
stat.meanTime=seq;
stat.fonum=seq;
rml=[];
for mi=1:markAmt
    % ��ǰmarker neuron spike���� marklb==1��ѡ����������Щ����ȫSTɨ��
    % * ST��������ҲҪ�������� ��ָֻ������Щmarklb==1��spike��
    sd=nsd2{mi}(marklb{mi}==1);
    [Ichi,Isd,rnum]=cspair(sd,ST,ext);
    
    % ����followerͨ������ֻ�����㹻�ظ�������ͨ��    
    if bAutoThres
        pthr=1e-6;
        
        foI=false(chAmt,1);
        for m=1:chAmt            
            [P,expectnum]=probable_rs(info.TimeSpan,[sAmt(seqMark(mi)),sAmt(m)],ext,rnum(m));
            % * when 2ch have spike 10000&1000, expect~30; when have
            % 10000&10000, expect~300-400.
%             if isnan(P)
%                 error('check this out: P calculation numeric problem');
%             end
            if rnum(m)>expectnum && P<pthr 
                % ��һ��������Ϊ�˱�֤P�㹻С��ԭ������Ϊrnumλ�ڸ��ʷֲ���ֵ����һ�ࣨ������󳬹����ģ���µ�Ԥ�ڣ���
                % �����෴������ԶС��Ԥ�ڣ���
                foI(m)=true;
            end            
        end
        foI=find(foI);
    else
        foI=find(rnum>=minThres);
    end
    
    foAmt=length(foI);
    Ichi=Ichi(:,foI); Isd=Isd(foI);
    if foAmt<=1, rml=[rml,mi]; end % ���û��follower,���ͨ��������seq,ȡ����
    % seq channels composition
    tp=(1:chAmt); seq{mi}=tp(foI);
    stat.fonum{mi}=rnum(foI); % number of followers in each channel
    
    % Average time
    mt=zeros(foAmt,1);
    for k=1:foAmt
        mt(k)=mean(ST{foI(k)}(Isd{k}));
    end
    % sort by time
    [mt,I]=sort(mt,'ascend');
    seq{mi}=seq{mi}(I);
    seqCS{mi}=Isd(I);
    stat.fonum{mi}=stat.fonum{mi}(I);
    stat.meanTime{mi}=mt;
end
if ~isempty(rml)
    seq(rml)=[];
    seqMark(rml)=[];
    seqCS(rml)=[];
    stat.markSDnum(rml)=[];
    stat.meanTime(rml)=[];
    stat.fonum(rml)=[];
    marklb(rml)=[];
end
stat.marklb=marklb;

end % of main