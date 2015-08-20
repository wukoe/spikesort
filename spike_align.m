%   [A,rmlist,O]=spike_align(X,SD,srate,varargin)
% X could be either signal data or the handle to .mat file.
% here SD is index form (not time format);
% Options:
% 'window',[pre,post]: specify the window size. (default: -0.8,1.1 ms)
% 'chAssign',chID: deal with when X and SD does not have one-to-one 
% correlations (like after sorting). chID is X-ch index of each SD channel.
% 'smooth',true/false: to smooth the edges of spikes.
% 'up sample',['up','down']: 
% Output:
% O:preww, postww, spklen, srate(in case up sampling)
function [A,rmlist,O]=spike_align(X,SD,srate,varargin)
%%% Default setting
% window for spike waveform
preww=0.8; % pre-event window width (ms)
postww=1.1; % post-event
% �Ƿ�ƽ���������������spike��ɵķ�
bSmoothOthers=false;
presmwin=[-1,-0.5]; % pre-peak smooth curve window - ����������ʼ�㵽�ź��յ�
% * ��ֵ��0.5s��1s��������1���ɵ�0.1
postsmwin=[0.5,1]; % post-peak - �ź���ʼ�㵽���������յ�
% Up sampling opt
usR=10;

%%% User input
bNewChID=false;
bUpSamp=false;
if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'smooth'
                bSmoothOthers=pinfo{parai};
            case 'window'
                preww=-pinfo{parai}(1);
                postww=pinfo{parai}(2);
            case 'chAssign'
                chID=pinfo{parai};                
                bNewChID=true;
            case 'up sample'
                bUpSamp=true;
                upsOpt=pinfo{parai};
            otherwise
                error('unidentified options');
        end
    end
end

%%% Proc
% Determine data type
if isa(X,'matlab.io.MatFile') % �����ļ�����    
    dataType=1;
    [pntAmt,XchAmt]=size(X,'X');
elseif isnumeric(X) % ֱ���������   
    dataType=2;
    [pntAmt,XchAmt]=size(X);
else 
    error('invalid data type');
end

% Check new channel ID assigned (when no assign info is provided, SD and X channel number must agree).
SDchAmt=length(SD);
if ~bNewChID 
    assert(SDchAmt==XchAmt,'SD and X channel number mismatch, and new chID is not assigned');
    chID=(1:SDchAmt)';
end
sAmt=cellstat(SD,'length');

% 
preww=round(preww/1000*srate);
postww=round(postww/1000*srate);
ptsAmt=preww+postww+1;
if bUpSamp
end

%%% Spike wave aligned to obtain window of data
A=cell(SDchAmt,1);
rmlist=zeros(0,2);
for chi=1:SDchAmt
    if sAmt(chi)>0
        if dataType==1
            chx=X.X(:,chID(chi));
        else
            chx=X(:,chID(chi));
        end
        
        % Discard the spike too close to the edge of recording that window does not fit
        % * only need to consider the two ends
        rlt=[];
        for k=1:min(3,sAmt(chi))  
            if SD{chi}(k)-preww<1
                rlt=[rlt; k];
            end
            if SD{chi}(end-k+1)+postww>pntAmt
                rlt=[rlt; sAmt(chi)-k+1];
            end
        end
        rlt=sort(rlt,'ascend');
        % Remove possible repeated index (in case too few spikes)
        if rlt>1
            I=find(diff(rlt)==0);
            rlt(I+1)=[];
        end        
        rla=length(rlt);
        if rla>0
            SD{chi}(rlt)=[];  sAmt(chi)=sAmt(chi)-rla;
            rmlist=[rmlist; ones(rla,1)*chi,rlt];
        end
        
        % Data segmentation (different methods)        
        A{chi}=zeros(ptsAmt,sAmt(chi));        
        for si=1:sAmt(chi)
            dslow=SD{chi}(si)-preww; dshigh=SD{chi}(si)+postww;
            
            % 1/N Naive.
            temp=chx(dslow:dshigh);
%             temp=temp-mean(temp);
%             temp=temp-X(I(si),chi); % make Y-axis align at peak too.
            % 2/N XY-align at mean-crossing point
%             temp=temp-mean(temp);
%             temp=(temp<0);
%             idx=find(temp(I(si)-dslow:end),1);
%             dslow=I(si)+idx-fix(ww/2);
%             dshigh=I(si)+idx+fix(ww/2);
%             temp=X(dslow:dshigh,chi);            
            % 3/N XY-align
%             vmax=X(I(si),chi);
%             [vmin,idx]=min(X(I(si):I(si)+120,chi));
%             idx=I(si)+round(idx/2);
%             temp=X(idx-ww:idx+ww,chi);
%             temp=temp-(vmax+vmin)/2;
            % e/N
            
            % Save align segment.
            A{chi}(:,si)=temp;
        end
    end
end

%%% Smoothing other peaks within the window
% * Method: 
% �Ϸ����� replace (part of)original signal with its sqrt root 
% �·������̶�ʱ����ڵı���ϵ��ƽ�����ɡ�
if bSmoothOthers
    %%% Prepare length parameters and smooth coef to be used    
    presmwin=fix(presmwin/1000*srate); 
    presmwinRawLen=presmwin(2)-presmwin(1)+1;
    postsmwin=fix(postsmwin/1000*srate); 
    postsmwinRawLen=postsmwin(2)-postsmwin(1)+1;
    
    %%% ͶӰ��[0,1] -> ��������ϵ�� % ��sigmoid
    tp=(0:presmwinRawLen-1)'/(presmwinRawLen-1);
    % 1/N
%     presmcoef=transition_curve(tp,'sigmoid','rise');
    % 2/N
    presmcoef=transition_curve(tp,'line','rise');
    % e/N
    presmcoef=presmcoef*0.9+0.1; % ��0-1������0.1-1    
    
    tp=(0:postsmwinRawLen-1)'/(postsmwinRawLen-1);
    % 1/N
%     postsmcoef=transition_curve(tp,'sigmoid','fall');    
    % 2/N
    postsmcoef=transition_curve(tp,'line','fall');
    % e/N
    postsmcoef=postsmcoef*0.9+0.1;
    
    % ת��Ϊpresmwin, postsmwin��A�����ж�Ӧ��λ�� 
    presmwin=presmwin+preww;
    postsmwin=postsmwin+preww;
    
    % ���Ȳ���ʱҲҪ������䲢�޸�presmwin & postsmwin��ֵ
    if presmwin(1)>1
        presmcoef=[0.1*ones(presmwin(1)-1,1); presmcoef];
        presmwin(1)=1;
    end
    if postsmwin(2)<ptsAmt
        postsmcoef=[postsmcoef; 0.1*ones(ptsAmt-postsmwin(2),1)];
        postsmwin(2)=ptsAmt;
    end
    
    % ��ֹ����<=0������
    if presmwin(1)<1
        presmcoef=presmcoef(1-presmwin(1)+1 : end);
        presmwin(1)=1;
    end
    % ��ֹ���ֳ���A���ݳ��ȵ�����
    if postsmwin(2)>ptsAmt
        postsmcoef=postsmcoef(1 : postsmwinRawLen-(postsmwin(2)-ptsAmt));
        postsmwin(2)=ptsAmt; 
    end    
    
    % 2/N ��0.1-1 change to 10-1
    postsmcoef=20-21*(postsmcoef-0.1);    
    presmcoef=20-21*(presmcoef-0.1);
    
    %%%
    preDlen=presmwin(2)-presmwin(1)+1; postDlen=postsmwin(2)-postsmwin(1)+1;
    for chi=1:seleAmt
        if sAmt(chi)>=1
            D=A{chi}(presmwin(1):presmwin(2),:);            
            for si=1:sAmt(chi)
                % 1/N
%                 D(:,si)=D(:,si).*presmcoef;
                % 2/N                
                tr=max(abs(D(:,si))); D(:,si)=D(:,si)/tr;
                for di=1:preDlen
                    D(di,si)=dpol(D(di,si),presmcoef(di));
                end
                D(:,si)=D(:,si)*tr;                
                % e/N
            end            
            A{chi}(presmwin(1):presmwin(2),:)=D;

            D=A{chi}(postsmwin(1):postsmwin(2),:);
            for si=1:sAmt(chi)
                % 1/N
%                 D(:,si)=D(:,si).*postsmcoef;
                % 2/N
                tr=max(abs(D(:,si))); D(:,si)=D(:,si)/tr;
                for di=1:postDlen
                    D(di,si)=dpol(D(di,si),postsmcoef(di));
                end
                D(:,si)=D(:,si)*tr;
                % e/N
            end
            A{chi}(postsmwin(1):postsmwin(2),:)=D;
        end
    end
end

%%% Up-sampling for making accurate alignment.
if bUpSamp
    UA=cell(SDchAmt,1); newA=cell(SDchAmt,1);
    for chi=1:SDchAmt
        [UA{chi},newA{chi}]=spkUpSamp(A{chi},'center',preww+1,'up fold',usR,'realign',[]);
    end
    if strcmp(upsOpt,'up') % up-sampled data
        A=UA;
    elseif strcmp(upsOpt,'down') % down sampling of the up-sampled (shifted)
        A=newA;
    else
        error('invalid option');
    end
end

%%% Output the length of the window
O=struct('spklen',ptsAmt);
if bUpSamp && strcmp(upsOpt,'up')
    O.srate=srate*usR;
    preww=preww*usR;
    postww=postww*usR;
end
O.preww=preww; O.postww=postww;
end %main


%%%%%%%%%%
% For depressing the outlier noise 
% when s=0.5, X=1 =roughly=>Y=1; when s=20, X=1 => Y=0.1;
function Y=dpol(X,s)
Y=1./(1+exp(-s*X))-0.5; 
Y=Y*4/s;
end