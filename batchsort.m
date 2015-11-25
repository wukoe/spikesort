% batch spike sort
% Go to the folder, and run: batchsort()
function bSuccess=batchsort(varargin)

% Whether to pass the file view input
if nargin>=1 && strcmp(varargin{1},'pass')
    bFileListWait=false;
else
    bFileListWait=true;
end

%%%%%%%%%%%% Search all qualified files in current folder.
%%% Find .mcd files
[fileList,fileAmt]=batchfile('*.mcd');
if fileAmt==0
    disp('no files found in current folder');
    bSuccess=true;
    return
end


%%% waiting for the "go" command
if bFileListWait
    seleI=(1:fileAmt);
    bNotDone=true;
    while bNotDone
        % Display file list
        fprintf('%d files found:\n',fileAmt);
        for fi=1:fileAmt
            fprintf('%d: %s\n',seleI(fi),fileList{seleI(fi)});
        end

        % wait for user
        s=input('Proceed? [y/n/s]:','s');
        if strcmp(s,'n')
            bSuccess=false; 
            fileList=[]; fileAmt=0;
            bNotDone=false;
        elseif strcmp(s,'y')
            bNotDone=false;
        elseif strcmp(s,'s') % select file.
            s=input('selection:','s');
            seleI=str2num(s);            
            fileAmt=length(seleI);
        end
    end
end
fileList=fileList(seleI);

%%%%%%%%%%%%%% Proceed the files found
% * in case sorting of any data is unsuccessful, using exception catching,
% and ignore those files.

% if matlabpool('size')==0
%     matlabpool open 6
% end

for fi=1:fileAmt    
    fprintf('Worker%d: Now #%d - %s\n',labindex,fi,fileList{fi});
%     try
    [~,fn,~]=fname(fileList{fi});
    spikesort(fn);%,'run step',{':','detect'});
        
%     catch ME % if failed
%         rethrow(ME);
%         infotext=sprintf('! %s has failed.\n',fileList{fi});
%         % print to both screen and log file
%         fprintf(1,infotext);
%         fprintf(fid,infotext);
%         
%         fprintf(1,ME.message);
%         fprintf(fid,ME.message);
%         fprintf('\n');
%     end
bSuccess=true;
end

% matlabpool close

disp('All files finished!');
