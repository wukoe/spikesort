%   [A,rmlist,O]=spike_align(X,SD,srate,varargin)
% X could be either signal data or the handle to .mat file.
% here SD is index form (not time format);
%   spike_align(...'chAssign',chID) deal with when X and SD does not have
% one-to-one correlations (like after sorting). chID is X-ch index of each
% SD channel.
%   spike_align(...'select',I) to using specified channel
% in this case, X keep the original data, SD should have the same number of
% channels as specified by I .
%   spike_align(...'bSmooth',true/false) to smooth the edges of spikes
%   ...'window' .. specify the window size. (default: -1 ~ 3 ms)
% O:preww, postww, spklen
function [A,rmlist,O]=spike_align(X,SD,srate,varargin)
%%% Default setting
% window for spike waveform
preww=0.8; % pre-event window width (ms)
postww=1.1; % post-event

% 是否平滑主峰两侧的其他spike造成的峰
bSmoothOthers=false;
presmwin=[-1,-0.5]; % pre-peak smooth curve window - 降噪作用起始点到信号终点
% 峰值后0.5s到1s，比例从1过渡到0.1
postsmwin=[0.5,1]; % post-peak - 信号起始点到降噪作用终点

%%% User input
bChSelect=false;
bNewChID=false;
if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'bSmooth'
                bSmoothOthers=pinfo{parai};
            case 'window'
                preww=-pinfo{parai}(1);
                postww=pinfo{parai}(2);
            case 'select'
                seleI=pinfo{parai};
                bChSelect=true;
            case 'chAssign'
                chID=pinfo{parai};                
                bNewChID=true;
            case 'T'
                T=pinfo{parai};
            otherwise
                error('unidentified options');
        end
    end
end

% Determine data type
if isa(X,'matlab.io.MatFile') % 输入文件类句柄    
    dataType=1;
    [pntAmt,XchAmt]=size(X,'X');
elseif isnumeric(X) % 直接输入变量   
    dataType=2;
    [pntAmt,XchAmt]=size(X);
else 
    error('invalid data type');
end

%%% Proc
% Check new channel ID assigned (when no assign info is provided, SD and X channel number must agree).
SDchAmt=length(SD);
if ~bNewChID
    if SDchAmt>XchAmt
        error('SD and X channel number mismatch, and new chID is not assigned');
    end
    chID=(1:SDchAmt)';
end

% Check the seleI number and contant.
seleAmt=length(SD);
if bChSelect
    if size(seleI,2)>size(seleI,1) % transpose the chI if necessary
        seleI=seleI';
    end
    
    I=(seleI<1) | (seleI>XchAmt);
    if any(I)
        error('invalid in chI: beyond the X range');
    end
    
    if length(seleI)~=seleAmt
        error('specified channel length not the same size as SD'); 
    end
end

% 
preww=round(preww/1000*srate);
postww=round(postww/1000*srate);

if ~bChSelect
    seleAmt=SDchAmt;
    seleI=1:seleAmt;    
end
A=cell(seleAmt,1);
sAmt=cellstat(SD,'length');
ptsAmt=preww+postww+1;


%%%%%%%%%%%%%%% Spike wave aligned to obtain window of data
rmlist=zeros(0,2);
for chCount=1:seleAmt
    chi=seleI(chCount);
    if sAmt(chCount)>0
        if dataType==1
            x=X.X(:,chID(chCount));
        else
            x=X(:,chID(chCount));
        end
        
        %%% Discard the spike too close to the edge of recording that window does not fit
        % * only need to consider the two ends
        rlt=[];
        for k=1:min(3,sAmt(chCount))  
            if SD{chCount}(k)-preww<1
                rlt=[rlt; k];
            end
            if SD{chCount}(end-k+1)+postww>pntAmt
                rlt=[rlt; sAmt(chCount)-k+1];
            end
        end
        rlt=sort(rlt,'ascend');
        % remove possible repeated index
        if rlt>1
            I=find(diff(rlt)==0);
            rlt(I+1)=[];
        end
        
        rla=length(rlt);
        if rla>0
            SD{chCount}(rlt)=[];  sAmt(chCount)=sAmt(chCount)-rla;
            rmlist=[rmlist; ones(rla,1)*chi,rlt];
        end
        
        %%%
        % Output init
        A{chCount}=zeros(ptsAmt,sAmt(chCount));        
        for si=1:sAmt(chCount)
            %%% data segmentation (different methods)
            dslow=SD{chCount}(si)-preww; dshigh=SD{chCount}(si)+postww;
            
%             % 1/N Y-align at maximum
%             temp=X(dslow:dshigh,chi);
%             temp=temp-X(I(si),chi);

            % 2/N Y-align no change
            temp=x(dslow:dshigh);
%             temp=temp-mean(temp);
            
%             % 3/N XY-align at mean-crossing point
%             temp=temp-mean(temp);
%             temp=(temp<0);
%             idx=find(temp(I(si)-dslow:end),1);
%             dslow=I(si)+idx-fix(ww/2);
%             dshigh=I(si)+idx+fix(ww/2);
%             temp=X(dslow:dshigh,chi);
            
%             % 4/N XY-align
%             vmax=X(I(si),chi);
%             [vmin,idx]=min(X(I(si):I(si)+120,chi));
%             idx=I(si)+round(idx/2);
%             temp=X(idx-ww:idx+ww,chi);
%             temp=temp-(vmax+vmin)/2;
            
            % save data in output            
            A{chCount}(:,si)=temp;
        end
    end
end


%%%%%%%%%%%% Smoothing other peaks within the window
% * Method: 
% 老方法： replace (part of)original signal with its sqrt root 
% 新方法：固定时间段内的比例系数平滑过渡。
if bSmoothOthers
    %%% Prepare length parameters and smooth coef to be used    
    presmwin=fix(presmwin/1000*srate); 
    presmwinRawLen=presmwin(2)-presmwin(1)+1;
    postsmwin=fix(postsmwin/1000*srate); 
    postsmwinRawLen=postsmwin(2)-postsmwin(1)+1;
    
    %%% 投影到[0,1] -> 过渡曲线系数 % 用sigmoid
    tp=(0:presmwinRawLen-1)'/(presmwinRawLen-1);
%     % 1/N
%     presmcoef=transition_curve(tp,'sigmoid','rise');
    % 2/N
    presmcoef=transition_curve(tp,'line','rise');
    % e/N
    presmcoef=presmcoef*0.9+0.1; % 从0-1提升到0.1-1    
    
    tp=(0:postsmwinRawLen-1)'/(postsmwinRawLen-1);
%     % 1/N
%     postsmcoef=transition_curve(tp,'sigmoid','fall');    
    % 2/N
    postsmcoef=transition_curve(tp,'line','fall');
    % e/N
    postsmcoef=postsmcoef*0.9+0.1;
    
    % 转换为presmwin, postsmwin在A数据中对应的位置 
    presmwin=presmwin+preww;
    postsmwin=postsmwin+preww;
    
    % 长度不足时也要进行填充并修改presmwin & postsmwin的值
    if presmwin(1)>1
        presmcoef=[0.1*ones(presmwin(1)-1,1); presmcoef];
        presmwin(1)=1;
    end
    if postsmwin(2)<ptsAmt
        postsmcoef=[postsmcoef; 0.1*ones(ptsAmt-postsmwin(2),1)];
        postsmwin(2)=ptsAmt;
    end
    
    % 防止出现<=0的索引
    if presmwin(1)<1
        presmcoef=presmcoef(1-presmwin(1)+1 : end);
        presmwin(1)=1;
    end
    % 防止出现超出A数据长度的索引
    if postsmwin(2)>ptsAmt
        postsmcoef=postsmcoef(1 : postsmwinRawLen-(postsmwin(2)-ptsAmt));
        postsmwin(2)=ptsAmt; 
    end    
    
    % 2/N 从0.1-1 change to 10-1
    postsmcoef=20-21*(postsmcoef-0.1);    
    presmcoef=20-21*(presmcoef-0.1);
    
    %%%
    preDlen=presmwin(2)-presmwin(1)+1; postDlen=postsmwin(2)-postsmwin(1)+1;
    for chCount=1:seleAmt
        if sAmt(chCount)>=1
            D=A{chCount}(presmwin(1):presmwin(2),:);            
            for si=1:sAmt(chCount)
%                 % 1/N
%                 D(:,si)=D(:,si).*presmcoef;
                % 2/N                
                tr=max(abs(D(:,si))); D(:,si)=D(:,si)/tr;
                for di=1:preDlen
                    D(di,si)=dpol(D(di,si),presmcoef(di));
                end
                D(:,si)=D(:,si)*tr;                
                % e/N
            end            
            A{chCount}(presmwin(1):presmwin(2),:)=D;

            D=A{chCount}(postsmwin(1):postsmwin(2),:);
            for si=1:sAmt(chCount)
%                 % 1/N
%                 D(:,si)=D(:,si).*postsmcoef;
                % 2/N
                tr=max(abs(D(:,si))); D(:,si)=D(:,si)/tr;
                for di=1:postDlen
                    D(di,si)=dpol(D(di,si),postsmcoef(di));
                end
                D(:,si)=D(:,si)*tr;
                % e/N
            end
            A{chCount}(postsmwin(1):postsmwin(2),:)=D;
        end
    end
end

%%%%%%%%%%% Output the length of the window
O=struct('preww',preww,'postww',postww,'spklen',ptsAmt);
end %main


%%%%%%%%%%
% For depressing the outlier noise 
% when s=0.5, X=1 =roughly=>Y=1; when s=20, X=1 => Y=0.1;
function Y=dpol(X,s)
Y=1./(1+exp(-s*X))-0.5; 
Y=Y*4/s;
end