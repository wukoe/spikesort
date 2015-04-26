% "run-through" for operating the file piece by piece.
%  input: original file; output: detected spike trains
%   mc2mat(filename,varoption)
% in varoption: "overwrite" cover existing .mat files; "pass" pass over
% these; "server" to run in server mode.
function success=mc2mat(fileName,varargin)
%%%%%%%%%%%% Basic setting 
%%% parameter default
blockSize=0.1; % read block size (ratio to the available memory)
bOverWrite=false;
bPass=false;
bServer=false;

% Output running state (fid=fopen(): to file; fid=1: to screen)
fid=1;

success=false; % the success return flag initialtion

%%% Handle the input
if ~isempty(varargin)
    for k=1:nargin-1
        switch varargin{k}
            case 'overwrite'
                bOverWrite=true;
            case 'pass'
                bPass=true;
            case 'server'
                % For special treatment of current server version of mat
                % lab that does not contain the 'matfile' part.
                bServer=true;
            case 'logf' % the FID of log file
                fid=varargin{k+1};
            otherwise
                % dismiss, because it may be the fid after 'logf'.
        end
    end
end


%%%%%%%%%%%%%%% Input and output 
%%% Locate the input file 
% setting: transformed file use the same name and same location as input file
[fp,fn,fex]=fname(fileName);

% if the input file extension is not .txt, warning
if strcmp(fex,'.mcd')
    intype=1;
elseif strcmp(fex,'.txt')
    intype=2;
else
    error('warning: input file extension not txt type\n');
end

% check if file with same name as planned output already exist
outputFileName=fn;
if exist([fp,fn,'.mat'],'file')
    if bPass
        disp('output file already exist, pass');
        return
    end
    if ~bOverWrite 
        % If not chosen overwrite file by default
        bNotDone=true;
        while bNotDone
            s=input('output file already exist, over write?: ','s');
            if strcmp(s,'y') % goes on                
                bNotDone=false;
            elseif strcmp(s,'n') % extend the name by (1)
                outputFileName=[fn,' (1)'];
                bNotDone=false;
            elseif strcmp(s,'r') % rename the file
                outputFileName=input('new name: ','s');
                bNotDone=false;
            else %the loop carries on
                disp('input again');
            end
        end
        
    end
end


%%% Open input file
if intype==1
    inputFileName=[fp,fn,'.mcd'];
    % Set library
    if bServer
        %%%%%%
    else
        nsresult = ns_SetLibrary('F:\Dropbox\mylibs\matlabTool\neuroshareMCDL\nsMCDlibrary64.dll');
        if nsresult==-1
            error('locate library error');
        end
    end
    
else
    inputFileName=[fp,fn,'.txt'];
    infile=fopen(inputFileName,'r');
    if infile==-1
        error('open .txt file error');
    end
end

%%% Set output mat file for storing
outputFileName=[fp,outputFileName,'.mat'];

if intype==1
    
else  
end

end

function convres=mcd2mat(outoutFileName,inputFileName)
[nsresult,fileInfo]=ns_GetFileInfo(hfile);
    %%%%
    
    [nsresult, infile] = ns_OpenFile(inputFileName);
    if nsresult==-1
        error('open .mcd file error');
    end
end

function convres=mcdtxt2mat(outoutFileName,inputFileName)
if bServer
    outfile=struct('X',[],'T',[],'srate',1);
else
    outfile=matfile(outputFileName,'Writable',true);
end

%%%%%%%%%%%%% Preprocessing of file information and system information 
%%% Obtain the basic information of file 
% (this part is file type specific, use corresponding treatment for each)

info=mcdtxtfile(infile,0);
%
tp=info.lineAmt*info.lineLen/1024^3;
fprintf('Total data size: %d lines (%2.2g GB)\n',info.lineAmt,tp);
fprintf(fid,'Total data size: %d lines (%2.2g GB)\n',info.lineAmt,tp);

% Get the sampling rate information from the time column of data.
fseek(infile,info.dStLoc,-1);
% read 5 lines of data
D=mcdtxtfile(infile,1,info,5);
temp=D.T;

% Check these 5 reading of sample time to see if their interval is stable
% (checking time data quality)
% But first get only the former 4 digits afte decimal point, ignore rest
temp=round(temp*1e4)/1e4; % * !! use only round() but not floor(), it will cause prolem - might be a bug of system.
sg=diff(temp); % intervals calculated by difference
if sum(abs(diff(sg)))>0 % difference of intervals 
    disp('warning: interval in time column not match');
    % remove the smallest and largest
    [~,idx]=max(sg); sg(idx)=[];
    [~,idx]=min(sg); sg(idx)=[];
end
sg=mean(sg); % Use mean as evaluation

% Transform interval (ms) to srate
srate=1000/sg;
% Round srate to ignore the 3 digits before decimal point (just for
% precaution)
srate=round(srate/1e3)*1e3;
fprintf('sampling rate: %d\n',srate);
fprintf(fid,'sampling rate: %d\n',srate);


%%% Determine the memory usage quantity 
% get system spare memory quantity
try
    [~,mem]=memory();
    mem=mem.PhysicalMemory.Available;
catch
    % If system mem info is unavailble, suppose it is 1GB on PC, 5GB on
    % server
    if bServer
        mem=10*1024^3; % won't be much when multiplied by the ratio
    else
        mem=1*1024^3;
    end
end
fprintf('memory to be used: %0.2g GB\n',mem*blockSize/1024^3);
fprintf(fid,'memory to be used: %0.2g GB\n',mem*blockSize/1024^3);

% turn the block size unit to (line)
blockSize=floor(mem*blockSize/info.lineLen);
fprintf('The processing block size: %d\n',blockSize);
fprintf(fid,'The processing block size: %d\n',blockSize);

%%% write important information to output
outfile.srate=srate;
temp=struct('chAmt',info.chAmt,'chID',info.chID,'ptsAmt',info.lineAmt);
outfile.info=temp;


%%%%%%%%%%%%%% the piece-reading part
% % to choose specific start point
% startLine=floor(startLoc/100 * lineAmt);
% if startLine+readLineNum>=lineAmt
%     readLineNum=lineAmt-startLine;
% end
% tp=startLine * lineLen + fhe;
% fseek(fid,tp,-1);
% fgetl(fid);

% calculate how many blocks are needed
ptsAmt=info.lineAmt;
[bSeg,bAmt]=cutseg(ptsAmt,blockSize);
fprintf('total number of blocks: %d\n',bAmt);
fprintf(fid,'total number of blocks: %d\n',bAmt);

% move to location of read start
fseek(infile,info.dStLoc,-1);
% initialize the output
outfile.X=zeros(1,info.chAmt,'double');
outfile.T=zeros(1,1,'double');

% read loop
if ~bServer, H=waitbar(0,num2str(ftell(infile))); end
for bi=1:bAmt
    D=mcdtxtfile(infile,1,info,blockSize);
    outfile.X(bSeg(bi,1):bSeg(bi,2),:)=D.X;
    outfile.T(bSeg(bi,1):bSeg(bi,2),1)=D.T;
    
    if bServer % save in each loop if it is a server version
        save(outputFileName,'-v7.3','-struct','outfile','X');
    end
    % visualization of progress
    fprintf('|');
    if ~bServer, waitbar(bi/bAmt,H,num2str(ftell(infile))); end
end
fprintf('\n');
if ~bServer, close(H); end

%%% Post-processing
fclose(infile);
% make a final save
if bServer
    save(outputFileName,'-v7.3','-struct','outfile');
%     save(outputFileName,'-v7.3','-struct','outfile','srate','X','T');
end
success=true;
end
