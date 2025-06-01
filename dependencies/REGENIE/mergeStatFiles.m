function mergeStatFiles(summFiles, outname, n)
if nargin < 3
    n = [];
end

if ispc
    cmdpref = "wsl /bin/bash -c ";
else
    cmdpref = "";
end

[bash_dir, bash_name] = fileparts(outname);


if ispc
    summFiles = string(summFiles);
    summFiles = arrayfun(@makeWSLpath, summFiles);
    outname = makeWSLpath(outname);
end

% if ispc
%     for i = 1:numel(summFiles)
%         dos2unix(summFiles(i));
%     end
% end

if ~isempty(n)
    % Remove first line (header) from summary stat files
    for kk = 2:numel(summFiles)
        tailcmd = "'"+string(summFiles{kk})+"'>'"+summFiles{kk}+".tmp'"+...
            " && mv '"+summFiles{kk}+".tmp' '"...
            +string(summFiles{kk})+"'"; 
        if ispc
            [~, ~] = system(cmdpref + '"tail -n ' + n + tailcmd + '"');
        else
            [~, ~] = system(cmdpref + 'tail -n ' + n + tailcmd);
        end
%         if fid % may fail with wsl when files are in another dir
%             tailcmd = '"' + string(summFiles{kk}) + '">"' + ...
%                 makeWSLpath(string(summFiles{kk}) + '.tmp', true) + ...
%                 '"' + ' && ' + cmdpref + 'mv "' + summFiles{kk} + ...
%                 '.tmp"' + ' "' + string(summFiles{kk}) ...
%                 + '"'; 
%             system(cmdpref + 'eval $"tail -n ' + n + tailcmd + '"');
%         end
    end
end

try
    summFiles(summFiles == outname) = [];
catch
    disp(summFiles)
    disp(outname)
end
mergedSummStat = string(join(summFiles, "' '"'));

runbash("cat '" + mergedSummStat + "' > '" + outname + "'", bash_name + "mergeStatFiles", "verbose", false, "dir", bash_dir);
% if ispc
%     [~, ~] = system(cmdpref + '"cat ''' + mergedSummStat + "' > '"' + outname + '''"');
% else
%     [~, ~] = system(cmdpref + "cat '" + mergedSummStat + "' > '" + outname + "'");
% end
% if fid
%     system(cmdpref + 'eval $"cat "' + mergedSummStat + '" > "' + ...
%         makeWSLpath(outname, true) + '""');
% end
end % END