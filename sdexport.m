% to export spike data (SD) to text file
%  sdexport(fileName, SD)
% the output write one row as a channel.
function sdexport(fileName, SD)
if iscell(SD)
    chAmt=length(SD);
else
    chAmt=size(SD,2);
end

fid=fopen(fileName,'w+t'); % use "w+t" so /n is possible with windows machine
if iscell(SD)
    for chi=1:chAmt
        fprintf(fid,'%d\t',SD{chi});
        fprintf(fid,'\n');
    end
else
    for chi=1:chAmt
        fprintf(fid,'%d\t',SD(:,chi));
        fprintf(fid,'\n');
    end
end

fclose(fid);