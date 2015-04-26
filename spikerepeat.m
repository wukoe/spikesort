% for detect the multi signals from a single source.
%   [stat,seqlist,seqCS,seqT]=spikerepeat(ST,TimeSpan(s))
% ! For ST, use only accurate time from exactST() !
function [stat,seqlist,seqCS,seqT]=spikerepeat(ST,TimeSpan,varargin)
% standard to find cluster sequence.
ext=0.000026;
totaltimelimit=0.002; % full length
minThres=4/60*TimeSpan; % this is for final, pre-scan should use lower bound 
% *(suggestion: for seq with length L, use minThres/L).
bNoFurther=true; % 1=just find all repeat pattern, no complex template matching.

if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'ext'
                ext=pinfo{parai};
            case 'length limit'
                totaltimelimit=pinfo{parai};
            case 'complex'
                bNoFurther=false;
            otherwise
                error('unidentified options');
        end
    end
end

% Proc
chAmt=length(ST);

% Output - statistics of the sequences.
stat=struct('TimeSpan',TimeSpan,'neibTimeLimit',ext*2,'totalTimeLimit',totaltimelimit);
stat.chspknum=cellstat(ST,'length'); stat.chfr=stat.chspknum/TimeSpan;


%%%%%%%%%%%% Get all tightly-pacted clusters
% Put spikes of all channels together
EV=[]; M=[]; % EV is interval of time occupied by every spike; M is channel and index information of each.
for k=1:chAmt
    temp=[ST{k}-ext,ST{k}+ext];
    EV=[EV;temp]; % time window of every spike events
    M=[M;ones(stat.chspknum(k),1)*k,(1:stat.chspknum(k))']; % [channel ID, spike ID of that channel] data
end
[~,eol,~,eolt]=findoverlap(EV); % eol is index of events (index of row of EV and M)
% *注意：此处不能保证每个eol内的列表都是按原来的时间顺序排的，需要进行整理:
% If too few raw clusters, stop and return.
stat.rawCluNum=length(eol);
if stat.rawCluNum<=2
    disp('only <2 cases');
    seqlist=[]; seqCS=[]; seqT=[];
    return
end

% 将在一个cluster内的spike进行时间排序。
for k=1:stat.rawCluNum
    temp=EV(eol{k},1);
    [~,I]=sort(temp,'ascend');
    eol{k}=eol{k}(I);
end

%%% Filtering by total time length, get number of qualified clusters (cluAmt)
a=eolt(:,2)-eolt(:,1)-2*ext;
I=a<totaltimelimit;
eol=eol(I); 
cluAmt=length(eol);
stat.cluNumAtimefilt=cluAmt;

%%% Find cluster channel index and time
% chl store each qualified clusters as a series of temporally ordered 
% index of channel and index of spike of that channel (in the form of rows of [ch ID, spike ID], same as M)
chl=cell(cluAmt,1);
T=cell(cluAmt,1); % time of spikes 
for k=1:cluAmt
    chl{k}=M(eol{k},:);
    
    len=size(chl{k},1);
    T{k}=zeros(len,1);
    for m=1:len
        T{k}(m)=ST{chl{k}(m,1)}(chl{k}(m,2));
    end
end


%%%%%%%%%%%% Find unique sequences (all cluster with same channel sequence combine into one seq type)
seqlist=cell(0);
% IDX is index of chl, which is 1to1 corresponded to all elements of seqlist;
IDX=cell(0);

seqlist{1}=chl{1}(:,1);
IDX{1}=1;
for k=2:cluAmt
    % Check whether current sequence is in identified seqlist
    seqAmt=length(seqlist); % length of list of identified seq types
    flagInList=false;
    for m=1:seqAmt
        % if exist, update it to correponding member.
        if isequal(chl{k}(:,1),seqlist{m})            
            IDX{m}=[IDX{m},k];
            flagInList=true;
            break
        end
    end
    % If not exist, add it to list as a new member.
    if ~flagInList
        seqlist{seqAmt+1}=chl{k}(:,1);
        IDX{seqAmt+1}=k;        
    end
end
seqlist=seqlist';

% Number of found unique seqs & how many repeats each seq have & seq length.
seqAmt=length(seqlist);
seqcount=cellstat(IDX,'length')';
seqlen=cellstat(seqlist,'length');
% record these infor.
stat.rsNumO1=seqAmt;
stat.seqcountO1=seqcount; % save the raw seq count

%%% 把chl,T等信息按unique sequence type分入各个sequence中。
seqCS=cell(seqAmt,1); seqT=seqCS;
for k=1:seqAmt
    temp=chl(IDX{k});
    seqCS{k}=zeros(seqcount(k),seqlen(k));
    for m=1:seqcount(k)
        seqCS{k}(m,:)=temp{m}(:,2)';
    end
    temp=T(IDX{k});
    seqT{k}=zeros(seqcount(k),seqlen(k));
    for m=1:seqcount(k)
        seqT{k}(m,:)=temp{m}';
    end
end

% % <<< export raw seq
% seqlist=cell(length(chl),1);
% for k=1:length(chl)
%     seqlist{k}=chl{k}(:,1);
% end
% TDol=[]; seqCS=[]; seqT=[];
% return


%%%%%%%%%%%% Filter seqs with low firing.
% First some statistics, 达到不同count（重复出现次数）标准的seq的数量
% stat.rsNumO2=sum(seqcount>=2); stat.rsNumO3=sum(seqcount>=3);
stat.rsNumOM1=sum(seqcount>=TimeSpan*1/60);  
stat.rsNumOM6=sum(seqcount>=TimeSpan*0.1);
stat.rsNumOM30=sum(seqcount>=TimeSpan*0.5);
stat.rsNumOM60=sum(seqcount>=TimeSpan);
stat.rsNumOM180=sum(seqcount>=TimeSpan*3);

%%% <<< pre- return
if bNoFurther
    % Filtering
    saveI=seqcount>minThres;
    
    % Update according to filtering
    seqlist=seqlist(saveI);
    seqcount=seqcount(saveI);
    seqlen=seqlen(saveI);
    seqT=seqT(saveI);
    seqCS=seqCS(saveI);
    
    stat.seqcount=seqcount;
    stat.seqlen=seqlen;
    return
end

%%% 1st stage Filtering
% Specify the over-threshold seqs (根据不同的seq len设定不同的阈值)
% * 按方案1，这里的标准定得高些，重新扫描的时候再把其他带回来。
saveI=false(seqAmt,1);
for k=2:4
    I=(seqlen==k);
    idx=(seqcount(I)>=minThres/k);
    temp=saveI(I); temp(idx)=true;
    saveI(I)=temp;
end
% 超过5的一律采用和5相同的标准
I=(seqlen>=5);
idx=(seqcount(I)>=minThres/5);
temp=saveI(I); temp(idx)=true;
saveI(I)=temp;

% % Save the below-threshold sequences，用于后面的二次扫描 ！！
% remseq=seqlist(~saveI); % raw sequence remaining
% remseqcount=seqcount(~saveI);
% remseqlen=seqlen(~saveI);
% remseqT=seqT(~saveI);
% remseqCS=seqCS(~saveI);
% remseqAmt=length(remseq);

% Update according to filtering
seqlist=seqlist(saveI);
seqcount=seqcount(saveI);
seqlen=seqlen(saveI);
seqT=seqT(saveI);
seqCS=seqCS(saveI);
seqAmt=length(seqlist);


%%% Calculate the average spike time difference and the ratio of "outlier repeats"
% whose time difference are far away to this average.
TDmean=cell(seqAmt,1); TDol=zeros(seqAmt,1);
for si=1:seqAmt
    temp=diff(seqT{si},[],2);    
    % Get index of outlier - if one interval is larger, then the whole seq
    % is outlier.
    ol=false(seqcount(si),1);
    for k=1:seqlen(si)-1
        ol=ol | logical(outlier_detect(temp(:,k)));
    end
    TDol(si)=sum(ol);
    TDmean{si}=mean(temp(~ol));
end


%%%%%%%%%%%%% Go back and re-scan the raw clusters
% * to Find whether there is similarity between sequences.
if false
addseqcount=zeros(seqAmt,1);
addseqCS=cell(seqAmt,1);
addseqT=cell(seqAmt,1);
for si=1:remseqAmt
    % 对每个备选seq，找到所有与之相近的正式seq。
    M=zeros(seqAmt,1);
    for k=1:seqAmt
        % Find longest matched segments (insertion of one bit is tolerated)
        [pl,bn]=longestmatch(seqlist{k},remseq{si});
        if pl>=min(length(seqlist{k}),length(remseq{si}))-4 && pl<=max(length(seqlist{k}),length(remseq{si}))+4
            if bn<=4, M(k)=pl;  end            
        end
    end
    
    %%% Check if the time interval of newly-found fall within constraint.
    % <<<< to be added
    
    % If there are multiple matching, choose ...
    [tp,idx]=max(M);
    if tp>0
        % 如果有多个匹配达到最优标准，找最短的那个
        if sum(M==tp)>1
            I=find(M==tp);
            [~,idx2]=min(seqlen(I));
            idx=I(idx2);
        end
        
        addseqcount(idx)=addseqcount(idx)+remseqcount(si);        
        % ！注意！新序列和原序列的序列（包括序列长度）可能不一样，所以CS，T信息的结合需要改变。
        % 逐位分析备选seq的通道量
        %<<< how about seqCS and seqT 尚未想好
        tempCS=zeros(remseqcount(si),seqlen(idx)); tempT=tempCS;
        for k=1:seqlen(idx)
            % search for its loc in seq.
            idx2=find(remseq{si}==seqlist{idx}(k),1);
            if isempty(idx2)
                tempCS(:,k)=NaN;
            else
                tempCS(:,k)=remseqCS{si}(:,idx2);
                tempT(:,k)=remseqT{si}(:,idx2);
            end
        end
        addseqCS{idx}=[addseqCS{idx};tempCS];
        addseqT{idx}=[addseqT{idx};tempT];
    end
end
stat.addseqcount=addseqcount;
stat.initseqcount=seqcount;
% Add new number to seq repeat number
seqcount=seqcount+addseqcount;
% <<< seqCS and seqT
end


%%% Now do the second stage of filtering
sAmt=length(seqlist);
saveI=false(sAmt,1);
% % 1/2 By number
% % Specify the over-threshold seqs
% for k=2:4
%     I=(seqlen==k);
%     idx=(seqcount(I)>=minThres/(k-1)); % 注意这里不一样
%     temp=saveI(I); temp(idx)=true;
%     saveI(I)=temp;
% end
% I=(seqlen>=5);
% idx=(seqcount(I)>=minThres/4);
% temp=saveI(I); temp(idx)=true;
% saveI(I)=temp;

% 2/2 by probability
% * After testing, this gives very tiny number of probability, thus not
% doing any effective filtering. (reason might be that the "uniform
% distribution" assumption dose not apply in real firing).
for k=1:sAmt
    tp=ceil(stat.chspknum(seqlist{k}));
    p=probable_rs(TimeSpan,tp,ext,seqcount(k));
    if p<1e-6
        saveI(k)=true;
    end
end

% e/2

% Remove Un-qualified
seqlist=seqlist(saveI);
seqcount=seqcount(saveI);
seqT=seqT(saveI);
seqCS=seqCS(saveI);
stat.seqcount=seqcount;
stat.seqlen=seqlen;

end