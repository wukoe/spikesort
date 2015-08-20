% Upscale the spikes sampling by spline. And re-align the peaks.
%   [Y,newA]=spkUpSamp(A,varargin)
% A= bunch of spikes from same source(unit) [len X spkAmt].
% 'center',x : specify the center of search region at x.
% 'up fold',r
% 'realign',[]: 
% 'realign',newcenteridx: 
function [Y,newA,stat]=spkUpSamp(A,varargin)
usR=10; % scale of up-sampling.
centR=1; % search real peak in range [cent-centR, cent+centR].

% User
flagSearchCenter=true;
bRealign=false;
if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'center'
                cent=pinfo{parai};
                flagSearchCenter=false;
            case 'up fold'
                usR=pinfo{parai};
            case 'realign'
                bRealign=true;
                if isempty(pinfo{parai})
                    bSearchNewloc=true;
                else
                    bSearchNewloc=false;
                    newcenteridx=pinfo{parai};
                end
            otherwise
                error('unidentified options');
        end
    end
end
if bRealign && ~bSearchNewloc
    flagSearchCenter=false;
end

% Proc
[spklen,sAmt]=size(A);
newspklen=(spklen-1)*usR+1;
% specify center of search by most frequent maximum value of A.
if bRealign
    if flagSearchCenter % by finding most frequent peak location.
        [~,I]=max(abs(A));
        lb=reabylb(I);
        [~,idx]=max(lb.typeAmt);
        cent=lb.types(idx);
    end
    % map cent and centR to up-sampled scale.
    if bSearchNewloc
        centR=centR*usR;
        cent=cent*usR+1; %*注意此处用于总长度非newspklen，而是有向两边延伸的binX。
    end
    % time axis for spline for up-sampling.
    binTime=0:spklen+1;
    xx=linspace(0,spklen+1,newspklen+2*usR)';
else
%     cent=(cent-1)*usR+1; 
    % time axis for spline for up-sampling.
    binTime=1:spklen;
    xx=linspace(1,spklen,newspklen)';
end

%%%
Y=zeros(newspklen,sAmt);
for si=1:sAmt
    binY=A(:,si);
    if bRealign
        binY=[0;binY;0]; % extend at both ends.
        % * extension is necessary, as new peak location will shift within [cent-centR, cent+centR].
        % pre-A extension.
        tp=binY(2)-binY(3);
        if tp*binY(2)>0
            binY(1)=binY(2)+tp/2;
        else
            binY(1)=binY(2)+tp;
        end
        % post-A extensionm.
        tp=binY(end-1)-binY(end-2);
        if tp*binY(end-1)>0
            binY(end)=binY(end-1)+tp/2;
        else
            binY(end)=binY(end-1)+tp;
        end
    end
    
    % Up sampling
    yy=spline(binTime,binY,xx);
    
    % Select right center.
    if bRealign
        if bSearchNewloc
            temp=yy(cent-centR:cent+centR);
            [~,newcenteridx]=max(abs(temp));
        end
        yy=yy(newcenteridx : newcenteridx+newspklen-1);
    end
    
    Y(:,si)=yy;
end

% Down-sampling the data to the same length as input (this time the maximum
% is at tip of peak.)
I=((1:spklen)'-1)*usR+1;
newA=Y(I,:);

% Output
stat.newSpkLen=newspklen;
if bRealign
    stat.newCentIdx=newcenteridx;
end

end