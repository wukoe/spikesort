% Look for the closely-intervaled spike pairs between specified channels
% ��Ե�ǰͨ������marksd������Ϊwd��segment��ʽ������������ͨ����follower�ĳ���λ�á�
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
            % ��һ��������Ϊ�˱�֤����һ����ʵ��P�㹻С��ԭ������Ϊrnumλ�ڸ��ʷֲ���ֵ���Ҳࣨ������󳬹����ģ���µ�Ԥ�ڣ���
            % �����෴������ԶС��Ԥ�ڣ���
            foI(chi)=true;
        end
    else
        if rpnum>=param.minThres
            foI(chi)=true;
        end
    end    
    % ��û������������pair number�����channel(chi)��������ص�CS������.
    if ~foI(chi)
        continue
    end
        
        
    %%%%%%%%%%%%%%%%% Filtering difference if required.
    % Ŀ�ģ�ֻ������ӽ�ƽ��ʱ���spike.
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
        % *������Ҫ����������������Ϊ�Ѿ�û���µ�foI(chi)�ᱻ��Ϊtrue,ֻ��ԭ��Ϊtrue�Ļ���false��
        % ����Ȼʹ��ԭ���Ļ�ʹ�ò��ϸ��chi��Ȼ����true.
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