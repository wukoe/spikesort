% Look for the closely-intervaled spike pairs between specified channels
% 相对当前通道（由marksd代表，化为wd的segment形式），其他所有通道的follower的出现位置。
%   [Iw,Isd,num]=cspair(marksd,ST,ext,markSDtotal)
%   [Iw,Isd,num,]=cspair(marksd,ST,ext,markSDtotal,param)
% markSDtotal is the number of all spikes which marksd is in.
function [IW,ISD,rpNum,varargout]=cspair(marksd,ST,ext,sAmtMark,param)
flagFilt=true;

% Proc
markAmt=length(marksd);
cha=length(ST);
sAmt=cellstat(ST,'length');
if markAmt==0
    IW=[]; ISD=cell(cha,1);
    return
end

%%%
wd=[marksd-ext,marksd+ext];% Make into bins
IW=false(markAmt,cha);
ISD=cell(cha,1);
foI=false(cha,1);% indicator of "follower" channel of marksd
for chi=1:cha
    % Find "follower".
    BI=binid(ST{chi},wd);    
    biAmt=cellstat(BI,'length');%how many bins each spike belongs to. length as ST{chi}.
    Isd=find(biAmt>0);
    Isdlen=length(Isd);

    % Number threshold.
    rpnum=Isdlen;
    if param.bAutoThres
        [P,expectnum]=probable_rs(param.timeSpan,[sAmtMark,sAmt(chi)],ext*2,rpnum);
        % * when 2ch have spike 10000&1000, expect~30; when have
        % 10000&10000, expect~300-400.
        if rpnum>expectnum && P<param.numPthr
            % 第一个条件是为了保证这样一个事实：P足够小的原因是因为rnum位于概率分布峰值的右侧（数量大大超过随机模型下的预期），
            % 而非相反（数量远小于预期）。
            foI(chi)=true;
        end
    else
        if rpnum>=param.minThres
            foI(chi)=true;
        end
    end    
    % 若没有满足条件的pair number，这个channel(chi)不存在相关的CS，放弃.
    if ~foI(chi)
        continue
    end
        
        
    %%%%%%%%%%%%%%%%% Filtering difference if required.
    % 目的：只保留最接近平均时差的spike.
    if flagFilt
        %%% Find the true mean Time difference.
        % All possible pairs of time difference.
        DT=[];
        for k=1:Isdlen
            tp=ST{chi}(Isd(k))-marksd(BI{Isd(k)});
            DT=[DT;tp];
        end
        % Locate the main region of most frequent difference by histogram.
        % use number in this bin to compute mean.
        [N,tp]=hist(DT);
        [~,idx]=max(N);
        histbinwidth=(tp(2)-tp(1))/2;% get (half) bin width.
        tprange=[tp(idx)-histbinwidth,tp(idx)+histbinwidth]; % most frequent's bin range.
        I = (DT>tprange(1)) & (DT<tprange(2));% identify samples within this range.        
        dtm=mean(DT(I));        
        % Recursively approaximate the "true mean", start from the center of best bin.
        maxround=2;
        for ri=1:maxround
            seleI=(abs(DT-dtm)<=param.DTJthr);
            dt=DT(seleI); 
            dtm=mean(dt); % update mean
        end 
        
        %%% Find qualified pairs by "true mean"
        for k=1:Isdlen
            temp=ST{chi}(Isd(k))-marksd(BI{Isd(k)});
            temp=abs(temp-dtm);
            I=find(temp<=param.DTJthr);
            
            % When spike belongs to multiple bins, leave one.
            if length(I)>1
                [~,idx]=min(temp(I));
                BI{Isd(k)}=BI{Isd(k)}(I(idx));
            elseif length(I)==1
                BI{Isd(k)}=BI{Isd(k)}(I);
            else
                BI{Isd(k)}=[];
            end            
        end
        % * now new biAmt should be {0,1} only.
        % Transform BI from cell to vector.
        temp=zeros(sAmt(chi),1);
        for k=1:sAmt(chi) %length(I)
            if ~isempty(BI{k})
                temp(k)=BI{k};
            end
        end
        BI=temp;
        
        % When there's multiple spikes in the same bin, leave one.
        lb=reabylb(BI);
        % ignore 0 type
        idx=find(lb.types==0);
        if ~isempty(idx)
            lb.types(idx)=[];
            lb.typeAmt(idx)=[];
            lb.ids(idx)=[];
            lb.cAmt=lb.cAmt-1;
        end        
        % 
        rpI=find(lb.typeAmt>1);
        if ~isempty(rpI) % means there's cases of multiple spikes in one bin.
            for k=1:length(rpI)
                temp=ST{chi}(lb.ids{rpI(k)})-marksd(BI(lb.types(rpI(k))));
                temp=abs(temp-dtm);
                [~,idx]=min(temp);
                
                I=lb.ids{rpI(k)};
                I(idx)=[];
                BI(I)=0;
            end
        end
        
        %%% Check for second time whether the remaining still fullfill number limit.
        % *这里需要把条件反过来，因为已经没有新的foI(chi)会被设为true,只有原来为true的会变成false。
        % 若依然使用原来的会使得不合格的chi依然保持true.
        rpnum=sum(BI>0);
        if param.bAutoThres
            [P,expectnum]=probable_rs(param.timeSpan,[sAmtMark,sAmt(chi)],ext*2,rpnum);
            if rpnum<expectnum || P>param.numPthr
                foI(chi)=false;
            end
        else
            if rpnum<param.minThres
                foI(chi)=false;
            end
        end
        
        %%% can also add shape consistency check here.<<<
    else
        % Transform BI from cell to vector when not filt done.
        temp=zeros(sAmt(chi),1);
        for k=1:sAmt(chi) %length(I)
            if ~isempty(BI{k})
                temp(k)=BI{k};
            end
        end
        BI=temp;
    end
    
    % Store data in Iw and ISD.
    I=(BI>0);
    ISD{chi}=I;
    IW(BI(I),chi)=true;
end

% Keep only effective CS channels's information.
IW=IW(:,foI); ISD=ISD(foI); 
rpNum=sum(IW);
foI=find(foI); 
varargout{1}=foI;

end % main