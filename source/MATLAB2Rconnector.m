function [checkElapsedTime, isFailed] = MATLAB2Rconnector(codeName, opts)
% MATLAB2Rconnector sends MATLAB jobs to R 
% 
% @13/09/2022: 'code' option was added for more convenient use.
% @17/12/2022: some bugs were fixed.
% @24APR2023: 'log' flag was added to save stderr and stdout to a log
% file.
% @26JUNE2023: 'wsl' falg was added to run the code in WSL.

arguments
    codeName {mustBeText}
    opts.delr (1,1) logical = true; % removes the r code along with Rout or RData files (if Rscript.exe is not used).
    opts.maxTime (1,1) double = inf; % max time allowed to be spent with the process in R (in sec)
    opts.Rhome {mustBeTextScalar}
    opts.printStdout (1,1) logical = false; % see processManager doc
    opts.code {mustBeVector} % codes to be written to codeName
    opts.log (1,1) logical = false
    opts.wsl (1,1) logical = false % to run the code in WSL (under Windows)
end

codeName = string(codeName);
codeName = regexprep(codeName, ".r$", "", "ignorecase");
codeName = codeName + ".r";

if isfield(opts, 'code')
    opts.code = string(opts.code);
    if isrow(opts.code), opts.code = opts.code'; end
    writematrix(opts.code, codeName, "FileType", "text", "QuoteStrings", "none")
end

if opts.wsl
    codeName_wsl = makeWSLpath(codeName);
    dos2unix(codeName_wsl)
    isFailed = runbash("nohup R CMD BATCH --no-restore --vanilla --slave " + codeName_wsl, ...
        regexprep(codeName, ".r$", ""));

else % Windows

    if ~isfolder(opts.Rhome) % wrong R version or directory
        % try to locate the R version
        files = {dir('C:\Program Files\R\').name};
        files = files(startsWith(files, 'R-'));
        files = natsort(files);
        if isempty(files)
            error("MATLAB2Rconnector cannot locate R.exe! Edit it manually!")
        elseif numel(files) > 1
            files = files{end};
        end
        opts.Rhome = 'C:\Program Files\R\' + string(files) + "\bin\x64";
    end
    
    if opts.printStdout || ~isinf(opts.maxTime)
        opts.Rhome = opts.Rhome + "\Rscript.exe";
    else
        opts.Rhome = opts.Rhome + "\R.exe";
    
    end

    opts.Rhome = string(opts.Rhome);
   
    
    if ~isinf(opts.maxTime) % there are serious issues with processManager when running in parallel
        proc = processManager('command', '"' + opts.Rhome + '" "' + codeName + '"', ...
            'printStdout', opts.printStdout, 'printStderr', opts.printStdout);
        tic
        while proc.running
            if ~isinf(opts.maxTime) && toc >= opts.maxTime% terminate the process if exceeds timeVal
                proc.stop;
                fprintf('MATLAB2Rconnector has exceeded maxTime for %s\n', codeName)
                pause(1)
            end
        end
        checkElapsedTime = toc;
    
        if proc.exitValue ~= 0
            disp("WARNING: R failed to run " + codeName)
            isFailed = true;
        else
            isFailed = false;
        end
        
    else
        tic
        if isscalar(codeName)
            if opts.printStdout
                [isFailed, ~] = system('"' + opts.Rhome + '" "' + codeName + '"', '-echo');
            else
                [isFailed, ~] = system('"' + opts.Rhome + '" CMD BATCH --vanilla --slave "' + codeName + '"');
            end
        else % run in the background
            batfile = regexprep(codeName, ".r$", ".bat");
            cmd = '"' + opts.Rhome + ...
                '" CMD BATCH --vanilla --slave "' + ...
                codeName + '"';
            cmd = ["("; 'start /B "" ' + cmd'; ") | pause"];
    
            writematrix(cmd, batfile(1), FileType="text", QuoteStrings=false)
            [isFailed, ~] = system("call " + batfile(1));
            delete(batfile(1))
        end
        checkElapsedTime = toc;
    end
end

codeName = string(codeName);
if opts.log
    for k = 1:numel(codeName)
        pth = fileparts(codeName(k));
        if isfile(codeName(k) + ".Rout")
            if isempty(pth) || pth == ""
                movefile(codeName(k) + ".Rout", fullfile(pwd, codeName + "log.txt"), "f")
            else
                movefile(codeName(k) + ".Rout", fullfile(pth, "log.txt"), "f")
            end
        end
    end
end

if opts.delr && ~isFailed
    
    for k = 1:numel(codeName)
        delete(codeName(k))
        if isfile(codeName(k) + ".Rout")
            delete(codeName(k) + ".Rout");
        end
        if isfile(".RData"); delete(".RData"); end
    end
end
end