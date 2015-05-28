% For make the "ST" variable from "_s.mat" file exported from spike
% detection process. This is a batch processing
%   spkdetect2st(fn) fn is usually '_s.mat','_c.mat'
%   (fn,thres) apply activity threshold to filter some channels (n/sec)
%   (...,'exact') force to use exact time.
%   (...,'new') overwrite existing files if any
function spkdetect2st(fn,varargin)
% Default
bThres=false;
bExact=false; % Whether to use exact time (peak by interpolation)
bKeepExist=true; % Append to existing _st.mat file or replace it

% User 
if nargin>1
    vn=length(varargin);
    for k=1:vn
        if isnumeric(varargin{k})
            bThres=true;
            actThres=varargin{k}; % (/s)
        elseif strcmp(varargin{k},'exact')
            bExact=true;
        elseif strcmp(varargin{k},'new')
            bKeepExist=false;
        end
    end
end

%%% Get qualified File list
if ischar(fn)
    % Search all qualified files in current folder if no specified file name
    [fileList,fileAmt]=batchfile(fn);
elseif iscell(fn)    
    fileList=fn;
    fileAmt=length(fileList);
else
    error('invalid file names');
end

fstr=fn(end-5:end);
flagS=false; flagC=false;
if strcmp(fstr,'_s.mat')
    flagS=true;
elseif strcmp(fstr,'_c.mat')
    flagC=true;
else
    error('unknow data names');
end

% Remove the post-fix 
tplen=6;
for fi=1:fileAmt
    fileList{fi}(end-tplen+1:end)=[];
end

%%% Do jobs
fprintf('%d files\n',fileAmt);
for fi=1:fileAmt
    fileName=fileList{fi};
    fprintf('%s\n',fileName);    
    
    load([fileName,'_f.mat'],'T');
    if bExact        
        load([fileName,'_f.mat'],'X');
    end
    if flagS
        load([fileName,'_s.mat'],'SD','info');
        SDTchID=info.chID;
        
        if bExact
            SDT=exactST(X,SD,T,info.srate,SDTchID);
        else
            SDT=idx2time(SD,T);
        end
        
        if bThres
            FTactThres=actThres*info.TimeSpan;
            sa=cellstat(SDT,'length');
            seleI=(sa>=FTactThres);
            SDT=SDT(seleI);
            SDTchID=SDTchID(seleI);
        end
        
        if exist([fileName,'_st.mat'],'file') && bKeepExist
            save([fileName,'_st.mat'],'-append','SDT','info','SDTchID');
        else
            save([fileName,'_st.mat'],'SDT','info','SDTchID');
        end
    end
    if flagC
        load([fileName,'_c.mat'],'NSD','info');
        STchID=info.chID;
         
        if bThres
            FTactThres=actThres*info.TimeSpan;
            sa=cellstat(NSD,'length');
            seleI=(sa>=FTactThres);
            NSD=NSD(seleI);
            STchID=STchID(seleI);
        end
       
        if bExact
            ST=exactST(X,NSD,T,info.srate,STchID);
        else
            ST=idx2time(NSD,T);
        end
        if exist([fileName,'_st.mat'],'file') && bKeepExist
            save([fileName,'_st.mat'],'-append','ST','info','STchID');
        else
            save([fileName,'_st.mat'],'ST','info','STchID');
        end
    end
end

end