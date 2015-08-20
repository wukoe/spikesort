% current polar figure of axon direction
%   cspathplot(opt,seqlist,stat,CL)
% opt='polar'
% opt='mea'
function varargout=cspathplot(opt,seqlist,stat,CL)
% Setting
repthr=100;
bUseSeqRepeat=false;

% Proc
seqAmt=length(seqlist);

%%%
if strcmp(opt,'polar')
    % Get polar representation of direction and length of CS. (only
    % measures direct line between head and tail) <<<< change this to
    % pieces later.
    eudis=csdirection(seqlist,CL);
    xd=eudis(:,1); yd=eudis(:,2);
    % if need to use polar 
    seqR=sqrt(xd.^2 + yd.^2);
    seqTheta=atan2(yd,xd);
    if bUseSeqRepeat
        seqR=seqR.*stat.seqcount;
    end
    
    %
    seqTheta=mod(seqTheta,pi*2);
    aspan=pi/8;
    DA=zeros(9,1);
    for di=1:7
        ar=[di*pi/4-aspan,di*pi/4+aspan];
        I=(seqTheta>=ar(1)) & (seqTheta<=ar(2));
        tp=seqR(I);
        if isempty(tp)
            DA(di+1)=0;
        else
            DA(di+1)=mean(tp);
        end
    end
    I=(seqTheta>=0) & (seqTheta<=aspan);
    I=I | ((seqTheta>=(2*pi-aspan)) & (seqTheta<=2*pi));
    tp=seqR(I);
    if isempty(tp)
        DA(1)=0;
    else
        DA(1)=mean(tp);
    end
    
    DA(9)=DA(1);
    ast=0:pi/4:2*pi; ast=ast';
    polar(ast,DA); % use "compass" format.
    
    %%%
elseif strcmp(opt,'mea')
    % Filt CS with enough repeats.
    for seqi=1:seqAmt
        rn=stat.fonum{seqi};
        I=rn>repthr;
        seqlist{seqi}=seqlist{seqi}(I);
    end
    % removal
    tp=cellstat(seqlist,'length');
    sele=find(tp>1);
    seqlist=seqlist(sele);
    fonum=stat.fonum(sele);
    seqAmt=length(seqlist); % new seqAmt
    
    % Get mean repeat number, path length(as number of channels involved).
    rep=zeros(seqAmt,1); plen=rep;
    for seqi=1:seqAmt
        rep(seqi)=mean(fonum{seqi});
        plen(seqi)=length(seqlist{seqi});
    end
    
    % Plot
    meapath(seqlist,CL);
    s=sprintf('num:%d,repeat:%d(%d),len:%f(%d)',seqAmt,mean(rep),max(rep),mean(plen),max(plen));
    title(s);
    varargout{1}=seqlist; varargout{2}=sele;
end

end