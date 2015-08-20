% find direction of the propagation of sequence
%   [HT,PL]=csdirection(seq,CL)
% Euclidean measure of X,Y distance of sequence (fold of grid)
% HT= head to tail,
% PL= each link in path.
function [HT,varargout]=csdirection(seq,CL,varargin)

% if ~isempty(varargin)
%     [pname,pinfo]=paramoption(varargin{:});
%     % process the parameter options one by one
%     for parai=1:length(pname)
%         switch pname{parai}
%             case ''
%             otherwise
%                 error('unidentified options');
%         end
%     end
% end
if isnumeric(seq)
    seq={seq};
end

seqa=length(seq);

% X,Y distance of head-tail of seq path. 
yd=zeros(seqa,1); xd=yd;
for seqi=1:seqa
    [r1,c1]=find(CL==seq{seqi}(1));
    [r2,c2]=find(CL==seq{seqi}(end));
    yd(seqi)=r2-r1; xd(seqi)=c2-c1;    
end
HT=[xd,yd];

% X,Y distance of each link in seq path. 
if nargout>1
    PL=cell(seqa,1);
    for seqi=1:seqa
        seqlen=length(seq{seqi});
        yd=zeros(seqlen-1,1); xd=yd;
        for k=1:seqlen-1
            [r1,c1]=find(CL==seq{seqi}(k));
            [r2,c2]=find(CL==seq{seqi}(k+1));
            yd(k)=r2-r1; xd(k)=c2-c1;
        end
        PL{seqi}=[xd,yd];
    end
    
    varargout{1}=PL;
end

end