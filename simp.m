% simple signal plot for a given period of signal
%   simp(X) or simp(X,T)
%   simp(matfile,startloc(s),length(s))
function varargout=simp(varargin)

%%%%%%%%%%% Input processing
% Check the input data type
if isnumeric(varargin{1})
    % ### if the input type is data matrix
    % - the input could be either just X or X with time T
    % Check for the time axis information 
    if nargin==1
        % When only X is available
        X=varargin{1};
        T=1:size(X,1);
    elseif nargin==2
        % When time info is added.
        X=varargin{1}; T=varargin{2};
        if size(T,2)>1
            error('time information should have only one column');
        end
    else
        error('invalid input term number');
    end
    
elseif strcmpi(class(varargin{1}),'matlab.io.MatFile')
    % ### if the input is matfile object
    % - need input of start location (by percentage), and length of
    % plotting epoch
    df=varargin{1};
    [ptsAmt,~]=size(df,'X');
    srate=df.srate;
    
    startp=varargin{2};
    if nargin==3
        readLen=varargin{3};
        if readLen>10 % maximum 10s of data
            disp('the data is too long for ploting!');
            return
        end
    else % show 1 second of data
        readLen=1;
    end
    startp=floor(startp*srate);
    if startp>ptsAmt
        disp('the start location is longer than the total length.');
        return
    end
    readLen=floor(readLen*srate);
    
    X=df.X(startp:startp+readLen,:);
    T=df.T(startp:startp+readLen,:);
    T=T/1000;
else
    error('unknown input data type');
end

%%%%%%%%%%% Plotting
Xstd=std(X);
dist=mean(Xstd)*8;
X=rmmean(X);
chAmt=size(X,2);

clf
hold on
for k=1:chAmt    
    plot(T,X(:,k)+dist*(chAmt-k),'b');
end
tp=T(end)-T(1);
axis([T(1)-0.05*tp,T(end)+0.05*tp,-dist*2,dist*(chAmt+2)]);
hold off

% Output
if nargout>0
    varargout{1}=X;
    varargout{2}=T;
end
