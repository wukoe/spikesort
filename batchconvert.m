% batch convert
% Go to the folder, and run: batchconvert()
function bSuccess=batchconvert(varargin)
%% Init
% determine the right symbol in the path according to OS
s=computer();
switch s
    case {'PCWIN64','PCWIN'}
        psymb1='\';
        bServer=false;
    otherwise % Linux
        psymb1='/';
        bServer=true;
end


% Whether to pass the file view input
if nargin>=1 && strcmp(varargin{1},'pass')
    bFileListWait=false;
else
    bFileListWait=true;
end


%% Search all qualified files in current folder.
fpath=pwd();
temp=dir();
% find those with .txt extension
fileAmt=length(temp);
count=0;
fileList=cell(1,1);
for fi=1:fileAmt
    if ~temp(fi).isdir
        [~,~,s]=fname(temp(fi).name);
        if strcmp(s,'.mcd')
            count=count+1;
            fileList{count}=temp(fi).name;
        end
    end
end
fileAmt=count;    

if fileAmt==0
    disp('no files found in current folder');
    bSuccess=[];
    return
end

% Display file list
fprintf('%d files found:\n',fileAmt);
% fprintf(logf,'%d files found:\n',fileAmt);
for fi=1:fileAmt
    fprintf('%s\n',fileList{fi});
end

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

%% Proceed the files found
for fi=1:fileAmt
    fprintf('Now #%d: %s\n',fi,fileList{fi});
%     fprintf(logf,'Now #%d: %s\n',fi,fileList{fi});
    try
        if bServer
            mc2mat([fpath,psymb1,fileList{fi}],'pass','server');
        else
            mc2mat([fpath,psymb1,fileList{fi}],'pass');
        end
    catch ME
        errtxt=ME.message;
        disp(errtxt);
%         fprintf(logf,errtxt);
%         rethrow(ME);
    end
end

bSuccess=true;
disp('All files finished!');