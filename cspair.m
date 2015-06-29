% Look for the closely-intervaled spike pairs between specified channels
% 相对当前通道（由marksd代表，化为wd的segment形式），其他所有通道的follower的出现位置。
%   [Iw,Isd,num]=cspair(marksd,ST,ext)
%   [Iw,Isd,num,]=cspair(marksd,ST,ext,param)
function [Iw,Isd,rpnum,varargout]=cspair(marksd,ST,ext,varargin)
marklen=length(marksd);
cha=length(ST);

if marklen==0
    Iw=[]; Isd=cell(cha,1); rpnum=zeros(1,cha);
    return
end

% Sometimes, >1 bins can overlap each other, leading to possible hosting
% the same spikes in different bins, delete this situation. 
D=diff(marksd);
tdI=find(D<ext*2);
for k=1:length(tdI)
    % Always delete the one that is closer to up/down stream spikes.
    % <<< do not consider 3 spikes together for now, only of pairs
    tp=[D(tdI(k)-1),D(tdI(k)+1)];
    if tp(1)>tp(2) % delete second spike
        marksd(tdI(k)+1)=-ext*2; % put it in location where no spike in ST can fall-in, instead of deleting directly.
    else
        marksd(tdI(k))=-ext*2;
    end
end

% Make into bins and find "follower".
wd=[marksd-ext,marksd+ext];
Iw=false(marklen,cha);
Isd=cell(cha,1);
for chi=1:cha
    [BI,BA]=binid(ST{chi},wd);
    Iw(:,chi)=(BA>0);
    Isd{chi}=(BI>0);
    
    % Also, Sometimes >1 ST{chi} spikes can be in the same bin, here decide to
    % delete all others except the most close one.
    eI=find(BA>1);
    if ~isempty(eI)
        for ei=1:length(eI)
            % marksd spike time of selected bin
            spktime=marksd(eI(ei));
            % find bin member in Isd
            I=find(BI==eI(ei));
            % time diff to the spike time
            dt=abs(ST{chi}(I)-spktime);
            % choose closest one
            [~,idx]=min(dt);
            I(idx)=[];
            Isd{chi}(I)=false;
        end
    end
end
rpnum=sum(Iw,1);


%%%%%%%%%%%%%%
if nargin>=4
    param=varargin{1};
    sAmtMark=length(marksd);
    sAmt=cellstat(ST,'length');
    
    %%% Number threshold.
    if param.bAutoThres
        sdnum=length(Isd);
        foI=false(sdnum,1);% indicator of "follower"
        for m=1:sdnum
            [P,expectnum]=probable_rs(param.timeSpan,[sAmtMark,sAmt(m)],ext*2,rpnum(m));
            % * when 2ch have spike 10000&1000, expect~30; when have
            % 10000&10000, expect~300-400.
            if rpnum(m)>expectnum && P<param.numPthr
                % 第一个条件是为了保证这样一个事实：P足够小的原因是因为rnum位于概率分布峰值的右侧（数量大大超过随机模型下的预期），
                % 而非相反（数量远小于预期）。
                foI(m)=true;
            end
        end
        foI=find(foI); foAmt=length(foI);
    else
        foI=find(rpnum>=param.minThres); foAmt=length(foI);
    end
    % Keep only effective CS channels's information.
    Iw=Iw(:,foI); Isd=Isd(foI); rpnum=rpnum(foI);
    
    % 若没有任何满足条件的conduction通道，这个channel(chi)不存在相关的CS，放弃.
    if foAmt==0
        varargout{1}=foI;
        return
    end
    
    %%% Filter spikes by Time difference consistence.
    for fi=1:foAmt
        %%% 1st round
        % Time difference between each pair.
        dt=marksd(Iw(:,fi))-ST{foI(fi)}(Isd{fi});
        % Get mean
        dtm=mean(dt);
        % Delete any pair whose difference deviate from this mean more than
        % threshold.
        I=(abs(dt-dtm)>param.DTJthr);
        tp=find(Iw(:,fi));
        Iwtemp=Iw(:,fi); Iwtemp(tp(I))=false;
        % * removing Isd must be placed later ?
        tp=find(Isd{fi});
        Isdtemp=Isd{fi}; Isdtemp(tp(I))=false;
        
        %%% 2nd time
        % Time difference between each pair.
        dt=marksd(Iwtemp)-ST{foI(fi)}(Isdtemp);
        % Get mean
        dtm=mean(dt);
        % Delete any pair whose difference deviate from this mean more than
        % threshold.
        I=(abs(dt-dtm)>param.DTJthr);
        tp=find(Iw(:,fi));
        Iw(tp(I),fi)=false;
        % * removing Isd must be placed later ?
        tp=find(Isd{fi});
        Isd{fi}(tp(I))=false;
        
        %%% can also add shape consistency check<<<
        
        % Update repeat number.
        rpnum(fi)=sum(Iw(:,fi));
    end
    
    % Check for second time whether the remaining still fullfill number limit.
    if param.bAutoThres
        seleI=false(foAmt,1);
        for fi=1:foAmt
            [P,expectnum]=probable_rs(param.timeSpan,[sAmtMark,sAmt(foI(fi))],ext*2,rpnum(fi));
            if rpnum(fi)>expectnum && P<param.numPthr
                seleI(fi)=true;
            end
        end
    else
        seleI=(rpnum>=minThres);
    end
    foI=foI(seleI);
    % Keep only effective CS channels's information.
    Iw=Iw(:,seleI); Isd=Isd(seleI); rpnum=rpnum(seleI);
    % 再次，若没有任何满足条件的conduction通道，放弃这个channel(chi).
    
    varargout{1}=foI;
end

end % main