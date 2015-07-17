% exact time of spike peaks.
%   [ST,SA]=exactST(X,SD,T,srate,chID)
%   [ST,SA]=exactST(A,SD,T,srate,mid)
% mid is center of seeking the time
function [ST,newSA]=exactST(X,SD,T,srate,varargin)
apSearchIntH=0.2; %(ms) -half of the interval to do spline

chAmt=length(SD);
sAmt=cellstat(SD,'length');
% flagX: whether input is raw signal X or aligned A.
if iscell(X)
    flagX=false;
    mid=varargin{1};
else
    flagX=true;
    if isempty(varargin)
        chID=1:chAmt;
    else
        chID=varargin{1};
    end
end

% Check channel numbers
if flagX
    assert(max(chID)<=size(X,2),'chID exceed channel number of X');
    assert(length(chID)==chAmt,'SD,chID cha not match');
else
    tp=cellstat(X,'size',2);
    assert(isequal(tp,sAmt),'X and SD not match');
    % tp=find(tp>0,1);
    % ALen=length(X{tp}(:,1));
end

% Proc
apSearchIntH=round(apSearchIntH/1000*srate);

%%%
ST=cell(chAmt,1); newSA=ST;
xpos=-apSearchIntH:apSearchIntH;
if flagX
    for chi=1:chAmt
        ST{chi}=zeros(sAmt(chi),1);
        newSA{chi}=zeros(sAmt(chi),1);
        for k=1:sAmt(chi)
            % 1/2 with X,T,SD
            try
                idxloc=SD{chi}(k)+xpos;
                binTime=T(idxloc);
                binX=X(idxloc,chID(chi));
                xx=linspace(binTime(1),binTime(end), apSearchIntH*20+1);
                yy=myspline(binTime,binX,xx);
            catch ME
                % if anything happens, mostly it is due to index out of
                % bound, then forget about the spike, move on to next.
                continue
            end
            
            if binX(apSearchIntH+1)>mean(binX)
                [tp,idx]=max(yy);
            else
                [tp,idx]=min(yy);
            end
            
            % assign new time and amplitude
            ST{chi}(k)=xx(idx);
            newSA{chi}(k)=tp;
        end
    end
else    
    for chi=1:chAmt
        ST{chi}=zeros(sAmt(chi),1);
        newSA{chi}=zeros(sAmt(chi),1);
        for k=1:sAmt(chi)
            temp=X{chi}(:,k);
            temp2=temp(mid+xpos);
            [~,idx]=max(abs(temp2));
            tp=mid+xpos;
            idx=tp(idx);
            binX=temp(idx+xpos);
            binTime=idx+xpos;
            xx=linspace(binTime(1),binTime(end), apSearchIntH*20+1);
            yy=myspline(binTime,binX,xx);
            
            if temp(idx)>mean(temp)
                [tp,idxx]=max(yy);
            else
                [tp,idxx]=min(yy);
            end
            
            % assign new time and amplitude
            dt=(xx(idxx)-idx)/srate;
            ST{chi}(k)=T(SD{chi}(k))+dt;
            newSA{chi}(k)=tp;
            
            %         %%% for visualization in debug
            %         cla
            %         hold on
            %         plot(binTime,binX);
            %         xp=xx(idxx);
            %         line([xp,xp],[tp-5e-9,tp+5e-9],'color','r');
        end
    end
end

end