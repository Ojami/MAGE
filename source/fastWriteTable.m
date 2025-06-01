function fastWriteTable(tab, opts)
% similar to write table but faster: uses a combination of tall/write (need
% parallel pool) and linux 'cat' and 'tail' commands. 
% Oveis Jamialahmadi, University of Gothenburg, May 2024.


arguments
    tab 
    opts.threads (1,1) {mustBeGreaterThanOrEqual(opts.threads, 2)} = 20 % otherwise use writetable
    opts.output {mustBeTextScalar} = fullfile(pwd, "table.txt") % output path and name
    opts.delimiter {mustBeTextScalar} = "\t"
end

assert(istable(tab) || strcmp(gather(classUnderlying(tab)), 'table'), "tab must be table or tall table!")
if isempty(gcp('nocreate'))
    parpool("Processes", opts.threads);
end


opts.output = string(opts.output);
% [~, ~, ext] = fileparts(opts.output);
% opts.output = regexprep(opts.output, ext + "$", ".txt");
opts.dir = fileparts(opts.output);
if opts.dir == ""
    opts.dir = pwd;
end

% temp dir for saving tall table
opts.tmp = fullfile(opts.dir, getRandomName("tmp", 6, sep="_"));
mkdir(opts.tmp);

matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.NullProgressReporter);
if ~istall(tab), tab = tall(tab); end
evalc('write(opts.tmp, tab, "FileType", "text", "Delimiter", opts.delimiter);');
matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.CommandWindowProgressReporter);

tmp_files = getfilenames(opts.tmp, "txt", "fullpath", true).txt;
tmp_files_wsl = arrayfun(@makeWSLpath, tmp_files);
cmd = "dos2unix " + tmp_files_wsl;
runbash(cmd, getRandomName("tmp", 6, sep="_"), "parallel", true);

cmd = string;
cmd(1, 1) = "cat " + tmp_files_wsl(1) + " > '" + makeWSLpath(opts.output) + "'";
cmd(2, 1) = "tail -q -n +2 " + join(tmp_files_wsl(2:end), " ") + " >> '" + makeWSLpath(opts.output) + "'";
runbash(cmd, getRandomName("tmp", 6, sep="_"), "parallel", false);
% system("wsl sudo cat " + join(tmp_files_wsl, " ") + " > tab3.txt");

% tmp_files = cellstr(tmp_files);
% delete(tmp_files{:});
rmdir(opts.tmp, "s")


%% UNDER DEVELOPMENT
% tic
% 
% cols = colnames(tab);
% dt = varfun(@(x)string(class(x)), tab);
% dt = dt{:, :};
% dtstr = ismember(dt.lower, ["categorical", "cell", "string"]);
% tab = convertvars(tab, dtstr, @string);
% % dt(dtstr) = "%s\r\n";
% % dt(~dtstr) = "%.15g\r\n";
% eolStr = "\n";
% 
% % write string table
% tmp_files = string; ct = 1;
% if any(dtstr)
%     cstr = cols(dtstr);
%     hdr_frmt = join(repmat("%s", numel(cstr), 1), "\t") + eolStr;
%     frmt = join(repmat("%s", numel(cstr), 1), "\t") + eolStr;
%     fid = fopen("tmp1.txt", 'w');
%     fprintf(fid, hdr_frmt, cstr);
%     fprintf(fid, frmt, tab{:, cstr}');
%     fclose(fid);
%     tmp_files(ct) = "tmp1.txt"; ct = ct + 1;
% end
% 
% if any(~dtstr)
%     cstr = cols(~dtstr);
%     hdr_frmt = join(repmat("%s", numel(cstr), 1), "\t") + eolStr;
%     frmt = join(repmat("%.15g", numel(cstr), 1), "\t") + eolStr;
%     fid = fopen("tmp2.txt", 'w');
%     fprintf(fid, hdr_frmt, cstr);
%     fprintf(fid, frmt, tab{:, cstr}');
%     fclose(fid);
%     tmp_files(ct) = "tmp2.txt"; 
% end
% 
% 
% % for k = 1:numel(cols)
% %     fid = fopen("tmp" + k + ".txt", 'w');
% %     % write header
% %     fprintf(fid,'%s\r\n', cols(k));
% %     fprintf(fid, dt(k), tab.(cols(k))');
% %     fclose(fid);
% % end
% % tmp_files = "tmp" + (1:numel(cols)) + ".txt";
% 
% 
% % for k = 1:numel(tmp_files), dos2unix(tmp_files(k)); end
% % for k = 1:numel(tmp_files)
% %     cmd(k, 1) = "tr -d '\015' < " + tmp_files(k) + " > " + tmp_files(k);
% % end
% system("wsl sudo paste -d '\t' " + join(tmp_files, " ") + " > tab3.txt");
% tmp_files = cellstr(tmp_files);
% delete(tmp_files{:});
% 
% toc

end % END