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


% ������������ڸĳ�ֱ����?mcd�ļ�����Ҫ�ˡ�
% %%% Ignore those proceed middle files
% % i.g., those with _f, _fl, _s, _a, _fea, _c, _fin  postfix.
% rmlist=[];
% for fi=1:fileAmt
%     [~,fn,~]=fname(fileList{fi});
%     if length(fn)>2
%         % find the last '_' character
%         idx=find(fn=='_');
%         if ~isempty(idx)
%             % all the letters from this and after
%             s=fn(idx(end):end);
%             temp=strcmp(s,{'_f';'_fl';'_s';'_a';'_c'});
%             if multilogic(temp,'or')
%                 rmlist=[rmlist,fi];
%             end
%         end
%     end
%     if length(fn)>4
%         temp=strcmp(fn(end-3:end),{'_fea';'_fin'});
%         if multilogic(temp,'or')
%             rmlist=[rmlist,fi];
%         end
%     end
% end
% 
% if ~isempty(rmlist)
%     fileList(rmlist)=[];
%     fileAmt=length(fileList);
% end

%%% Display file list
fprintf('%d files found:\n',fileAmt);
for fi=1:fileAmt
    fprintf('%s\n',fileList{fi});
end

% waiting for the "go" command
bNotDone=bFileListWait;
while bNotDone
    s=input('Proceed? [y/n]:','s');
    if strcmp(s,'n')
        bSuccess=[]; 
        return
    elseif strcmp(s,'y')
        bNotDone=false;
    % else loop again
    end
end


%%%%%%%%%%%%%% Proceed the files found
% * in case sorting of any data is unsuccessful, using exception catching,
% and ignore those files.

% if matlabpool('size')==0
%     matlabpool open 7
% end

for fi=1:fileAmt    
    fprintf('Worker%d: Now #%d - %s\n',labindex,fi,fileList{fi});
%     try
    [~,fn,~]=fname(fileList{fi});
    spikesort(fn,'RunStep',{'spike detect',':'});
        
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
end

% matlabpool close
bSuccess=true;
disp('All files finished!');
