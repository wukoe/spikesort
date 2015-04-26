% find direction of the propagation of sequence
%   [r,theta]=csdirection(seq,CL)
% Give the polar axis measure of direction.
%   [x,y]=csdirection(seq,CL,'measure','eu') 
% Change output to Euclidean measure.
%   (...'seq repeat', SR) multiply the length by the repeat number SR.
function varargout=csdirection(seq,CL,varargin)
opt='polar';
bUseSeqRepeat=false;

if ~isempty(varargin)
    [pname,pinfo]=paramoption(varargin{:});
    % process the parameter options one by one
    for parai=1:length(pname)
        switch pname{parai}
            case 'seq repeat'
                bUseSeqRepeat=true;
                seqcount=pinfo{parai};
            case 'measure'
                opt=pinfo{parai};
            otherwise
                error('unidentified options');
        end
    end
end

%%% X,Y distance of Start/End of seq path. 
sa=length(seq);
yd=zeros(sa,1); xd=yd;
for si=1:sa
    [r1,c1]=find(CL==seq{si}(1));
    [r2,c2]=find(CL==seq{si}(end));
    yd(si)=r2-r1; xd(si)=c2-c1;    
end
if bUseSeqRepeat
    yd=yd.*seqcount; xd=xd.*seqcount;
end

%%% Output format
if strcmp(opt,'eu')
    varargout{1}=xd; varargout{2}=yd;
elseif strcmp(opt,'polar')
    varargout{1}=sqrt(xd.^2 + yd.^2);
    varargout{2}=atan2(yd,xd);
else
    error('invalid opt');
end

end