function gcta(gctab, bedwd, opts)
% performs GCTA-COJO analysis. 
% example: gcta('E:\REGENIE\ALT\ALTn.BMI.interaction\ALTn_INT_SNPxBMI.txt', fullfile(pwd, 'imputed/'))

arguments
    gctab {mustBeFile}
    bedwd {mustBeFolder} % bed files directory
    opts.mafile
    opts.lof_flag (1,1) logical = false
    opts.p (1,1) double {mustBeInRange(opts.p, 0, 1)} = 5e-8
    opts.win (1,1) double {mustBePositive} = 1000 % in kbp 
    opts.collinear (1,1) double = 0.9 
    opts.parallel (1,1) logical = true % for readGWASfile only

    % this flag is used when parallel is true, but gcta jobs should NOT be
    % run within a parfor loop. In this manner, both 'threads' and
    % 'workers' can be set to reasonably large values, and at the same time
    % does not compromise the efficiency. 
    opts.gctaParfor (1,1) logical = false 
    opts.workers (1,1) double = 30 % for readGWASfile only
    opts.threads (1,1) double = 15 % for GCTA
    opts.gctawd {mustBeFolder} = fileparts(which("gcta.m")) % gcta binary software directory
    opts.output {mustBeTextScalar} % output file name, if left empty, will be same as the input 'gctab'
end

if ~opts.parallel
    opts.gctaParfor = false; % works only when 'parallel' is true
end

% check parallel pool
if opts.parallel && isempty(gcp('nocreate'))
    if isnan(opts.workers)
        % number of physical cores
        maxNumCompThreads('automatic'); % reset (only if was changed in current session)
        opts.workers = maxNumCompThreads; 
    end
    parc = parcluster('Processes');
    parc.NumWorkers = opts.workers;
    parpool(parc); % should I delete in the end? 
end

% locate bed files
genofiles = getfilenames(bedwd, "bed").bed;
opts.bgen = false;
if isempty(genofiles)
    % check for bgen files
    genofiles = getfilenames(bedwd, "bgen").bgen;
    if isempty(genofiles)
        error('gcta:cannot find any bed/bgen file within %s', bedwd)
    else
        % as of 30OCT2024: COJO still does not support bgen
        error("gcta:COJO still does not support bgen inputs!")
        opts.bgen = true;
    end
elseif numel(genofiles) < 22
    warning('gcta:found less than 22 bed files!')
end


if opts.bgen
    genofiles = regexprep(genofiles, '.bgen$', ''); 
else
    genofiles = regexprep(genofiles, '.bed$', '');
end
genoPattern = strshare(genofiles, pattern=true);

% read summary stat file
if isfield(opts, 'output')
    output = opts.output;
else
    [hm, output] = fileparts(gctab);
    output = string(fullfile(hm, output));
end
opts.outdir = string(fileparts(output));
if opts.outdir == "", opts.outdir = pwd; end
if ~endsWith(opts.outdir, filesep), opts.outdir = opts.outdir + filesep; end

gctab = readGWASfile(gctab, 'parallel', opts.parallel, ...
    'full', true, 'light', true, 'n', true);

% Since p-value cannot go under realmin, modify 0 p-values based on
% the corresponding Z-scores 
f_zero = gctab.p <= realmin;
if any(f_zero)
    Z = abs(gctab.se(f_zero)./gctab.beta(f_zero)); % 1/Z
    gctab.p(f_zero) = realmin.*(Z./max(Z));
end

gctab.pos = []; % not needed 

chru = unique(gctab.chr); % must be numeric (i.e. X chr should be 23)
name = getRandomName("cojo"); % temp GCTA files
genofiles = compose(genoPattern, double(chru));

fprintf('splitting the GWAS data into chromosome tables...')
t1 = tic;
chrtab = cell(numel(chru), 1);
for k = 1:numel(chru)
    idx = gctab.chr == chru(k);
    chrtab{k} = gctab(idx, :);
    gctab(idx, :) = [];
    
    if ~opts.bgen
        % Change SNP names based on varid in BIM files
        chrtab{k} = snptag(chrtab{k}, fullfile(bedwd, compose(genoPattern, double(chru(k)))), false);
    end
    
    % change column names
    chrtab{k}.chr = []; 
    cols = chrtab{k}.Properties.VariableNames;
    cols(ismember(cols, 'allele1')) = {'A1'};
    cols(ismember(cols, 'allele0')) = {'A2'};
    cols(ismember(cols, 'afreq')) = {'freq'};
    cols(ismember(cols, 'beta')) = {'b'};
    cols(ismember(cols, 'annotate')) = {'SNP'};
    cols(ismember(lower(cols), 'n')) = {'N'};
    chrtab{k}.Properties.VariableNames = cols;
    cols = {'SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N'};
    [~, idx] = ismember(cols, chrtab{k}.Properties.VariableNames);
    chrtab{k} = chrtab{k}(:, idx);

end
fprintf('\b\b done (%.0f sec)\n', toc(t1))

fprintf('writing to temp .ma files...')
t1 = tic;
outdir = opts.outdir;
if opts.parallel
    parfor k = 1:numel(chru)
        writetable(chrtab{k}, fullfile(outdir, name + "." + chru(k) + ".ma"), FileType='text', Delimiter=char(9));
        dos2unix(fullfile(outdir, name + "." + chru(k) + ".ma"))
    end
else
    for k = 1:numel(chru)
        writetable(chrtab{k}, fullfile(outdir, name + "." + chru(k) + ".ma"), FileType='text', Delimiter=char(9));
        dos2unix(fullfile(outdir, name + "." + chru(k) + ".ma"))
    end
end

fprintf('\b\b done (%.0f sec)\n', toc(t1))

opts.gctawd = makeWSLpath(opts.gctawd);
bedwd = makeWSLpath(bedwd);
if ~endsWith(bedwd, "/"); bedwd  = bedwd + "/"; end

genoline = strings(numel(genofiles), 1);
for k = 1:numel(genofiles)
    if opts.bgen
        genoline(k) = "--bgen " + bedwd + compose(genoPattern, chru(k)) + ".bgen --sample " + bedwd + compose(genoPattern, chru(k)) + ".sample";
    else
        genoline(k) = "--bfile " + bedwd + compose(genoPattern, chru(k));
    end
end

t1 = tic;
if opts.gctaParfor
    gctawd = opts.gctawd;
    p = opts.p;
    win = opts.win;
    threads = opts.threads;
    collinear = opts.collinear;
    outdir = opts.outdir;
    outdirw = makeWSLpath(outdir);

    parfor i = 1:numel(genofiles)
        fprintf('%d (of %d)-Running cojo for %s\n', i, numel(genofiles), compose(genoPattern, chru(i)))
        
        gcta_command = "wsl " + gctawd + "/gcta64 " + ...
            genoline(i) +...
            " --cojo-file " + outdirw + name + "." + chru(i) + ".ma" + ...
            " --cojo-slct --cojo-p " + p + ...
            " --cojo-wind " + win + " --cojo-collinear " + collinear + ...
            " --out " + outdirw + name + i + " --thread-num " + threads + ...
            " --cojo-actual-geno"; % --cojo-actual-geno seems to be deprecated in new GCTA 
    
        [~, ~] = system(gcta_command);
    
        if isfile(fullfile(outdir, name + i + ".jma.cojo"))
            cojo = readtable(fullfile(outdir, name + i + ".jma.cojo"), 'FileType', 'text', 'TextType', 'string');
            saveME(output + i + ".cojo.mat", cojo)
            delete(fullfile(outdir, name + i + ".jma.cojo"))
            delete(fullfile(outdir, name + i + ".cma.cojo"))
            delete(fullfile(outdir, name + i + ".ldr.cojo"))
        end
        delete(fullfile(outdir, name + "." + chru(i) + ".ma"))
        fprintf('job %d is done\n', i)
    end
else
    for i = 1:numel(genofiles)
        fprintf('%d (of %d)-Running cojo for %s\n', i, numel(genofiles), compose(genoPattern, chru(i)))

        gcta_command = "wsl " + opts.gctawd + "/gcta64 " + ...
            genoline(i) +...
            " --cojo-file " + makeWSLpath(opts.outdir) + name + "." + chru(i) + ".ma" + ...
            " --cojo-slct --cojo-p " + opts.p + ...
            " --cojo-wind " + opts.win + " --cojo-collinear " + opts.collinear + ...
            " --out " + makeWSLpath(opts.outdir) + name + i + " --thread-num " + opts.threads + ...
            " --cojo-actual-geno"; % --cojo-actual-geno seems to be deprecated in new GCTA 
    
        [~, ~] = system(gcta_command);
    
        if isfile(fullfile(opts.outdir, name + i + ".jma.cojo"))
            cojo = readtable(fullfile(opts.outdir, name + i + ".jma.cojo"), 'FileType', 'text', 'TextType', 'string');
            saveME(output + i + ".cojo.mat", cojo)
            delete(fullfile(opts.outdir, name + i + ".jma.cojo"))
            delete(fullfile(opts.outdir, name + i + ".cma.cojo"))
            delete(fullfile(opts.outdir, name + i + ".ldr.cojo"))
        end
        delete(fullfile(opts.outdir, name + "." + chru(i) + ".ma"))
        fprintf('job %d is done\n', i)
    end
end
fprintf('\nall jobs done (%.0f sec)\n', toc(t1))

% merge log files
logfis = getfilenames(opts.outdir, "log", "fullpath", false).log;
logfis(~startsWith(logfis, name)) = [];
logfis = fullfile(opts.outdir, logfis);
logout = readall(fileDatastore(logfis, 'ReadFcn', @readlines));
logout = vertcat(logout{:});
writelines(logout, output + ".cojo.log")
arrayfun(@delete, logfis)

% merge all mat files
matfiles = cellstr(output + (1:numel(genofiles)) + ".cojo.mat");
matfiles = matfiles(logical(cellfun(@exist, matfiles)));
cojo = readall(fileDatastore(matfiles, 'ReadFcn', @load));
cojo = vertcat(cojo{:});
cojo = struct2table(cojo, 'AsArray', true);
try
    cojo = vertcat(cojo.cojo{:});
catch
    cojo = cojo.cojo;
end

% there is a bug in 'gcta' wrapper function, where it doesn't convert snp
% names back to their original names
idx = startsWith(cojo.SNP, "rs") & contains(cojo.SNP, "_");
cojo.SNP(idx) = extractBefore(cojo.SNP(idx), "_");

save(output + ".cojo.mat", 'cojo')
delete(matfiles{:})

end % END 

%% subfunctions ===========================================================
function saveME(name, cojo)
    save(name, 'cojo')
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function makegcta(phenofile, bedRange, imp_flag, lof_flag, qvalue)
% 
% warning('imputed variants will have a modified varid in the output (snptag), consider changing it!')
% 
% if nargin < 5
%     qvalue = [];
% else
%     qvalue = string(qvalue);
% end
% if nargin < 4
%     lof_flag = false;
% end
% if nargin < 3
%     imp_flag = [];
% end
% if nargin < 2 || isempty(bedRange)
%     % chr identifier for imputed data
%     fprintf('ERROR:bedRange cannot be empty!\n')
%     return
% elseif numel(bedRange) > 1
%     fprintf('ERROR:bedRange length cannot be >1!\n')
%     return
% end
% 
% [~, rawName] = fileparts(phenofile);
% phenofile = load(phenofile);
% fname = fieldnames(phenofile);
% phenofile = phenofile.(fname{1});
% 
% if ismember('N', phenofile.Properties.VariableNames)
%     boltflag = false;
%     fprintf('SAIGE GWAS summary is provided.\n')
% else
%     boltflag = true;
%     fprintf('BOLT-LMM GWAS summary is provided.\n')
% end
% 
% gctab = table;
% 
% if boltflag
%     N = 404465; % Nused in BOLT-LMM log file (ALTn)
%     fprintf('Sample size: %d\n', N)
%     gctab.SNP = phenofile.SNP;
%     gctab.A1 = phenofile.ALLELE1;
%     gctab.A2 = phenofile.ALLELE0;
%     gctab.freq = phenofile.A1FREQ;
%     f_flip = gctab.freq > 0.5;
%     gctab.b = phenofile.BETA;
%     gctab.se = phenofile.SE;
%     try
%         gctab.p = phenofile.P_BOLT_LMM;
%     catch
%         gctab.p = phenofile.P_BOLT_LMM_INF;
%     end
%     try % imputed data do not contain this column
%         gctab.N = round(N.*(1-phenofile.F_MISS));
%     catch
%         gctab.N = N .* ones(size(gctab, 1), 1);
%     end
% else % SAIGE ----------------------
%     try
%         gctab.SNP = phenofile.rsid;
%     catch
%         gctab.SNP = phenofile.SNPID;
%     end
%     gctab.A1 = phenofile.Allele2;
%     gctab.A2 = phenofile.Allele1;
%     gctab.freq = phenofile.AF_Allele2;
%     gctab.b = phenofile.BETA;
%     if iscellstr(gctab.b)
%         gctab.b = double(string(gctab.b));
%     end
%     gctab.se = phenofile.SE;
%     try
%         gctab.p = phenofile.("p.value");
%     catch
%         gctab.p = phenofile.p_value;
%     end
%     gctab.N = phenofile.N;
%     if iscellstr(gctab.p)
%         gctab.p = double(string(gctab.p));
%     end
%     gctab(isnan(gctab.p), :) = [];
%     
%     f_flip = gctab.freq > 0.5;
% end
% 
% if any(f_fl ...
%         ip)
%     newA1 = gctab.A2(f_flip);
%     gctab.A2(f_flip) = gctab.A1(f_flip);
%     gctab.A1(f_flip) = newA1;
%     gctab.b(f_flip) = -1.*gctab.b(f_flip); % Binary traits must have log(OR)
%     gctab.freq(f_flip) = 1 - gctab.freq(f_flip);
% end
% 
% % Since p-value cannot go under realmin, modify 0 p-values based on
% % the corresponding Z-scores 
% f_zero = gctab.p <= realmin;
% if any(f_zero)
%     Z = abs(gctab.se(f_zero)./gctab.b(f_zero)); % 1/Z
%     gctab.p(f_zero) = realmin.*(Z./max(Z));
% end
% 
% % \imputed directory
% if isempty(imp_flag)
%     imp_flag = input('Imputed variants?[Y/n] ');
%     if strcmpi(imp_flag, 'n')
%         imp_flag = false;
%     else
%         imp_flag = true;
%         lof_flag = input('Only pLOF?[y/N] ');
%         if strcmpi(lof_flag, 'y')
%             lof_flag = true;
%         else
%             lof_flag = false;
%         end
%     end
% end
% 
% 
% if imp_flag
%     if lof_flag
%         imp_dir = 'imputed.lof';
%     else
%         imp_dir = 'imputed';
%     end
% 
%     bedfile = dir(fullfile(pwd, imp_dir, '*.bed'));
%     bedfile = string({bedfile.name});
%     bedfile = regexprep(bedfile, '.bed$', '');
%     bedfile = sortFiles(bedfile);
%     if ~lof_flag
%         if numel(bedfile) > 22
%             fprintf('ERROR:number of chromosomes canno be > 22!\n')
%             return
%         elseif numel(bedfile) < 22
%             numtag = double(string(regexp(bedfile, '\d+$', 'match')));
%             bedfile = bedfile(bedRange == numtag);
%         else % == 22
%             bedfile = bedfile(bedRange);
%         end
%     end
%   
%     % Change SNP names based on varid in BIM files
%     gctab = snptag(gctab, fullfile(imp_dir, bedfile), unique(phenofile.CHR));
% end
% clear phenofile
% 
% writetable(gctab, rawName+".pheno.ma", 'FileType', 'text', 'Delimiter', ' ');
% if isempty(qvalue)
%     qvalue = string(0.05/size(gctab, 1));
%     fprintf('qvalue: %s',qvalue)
% 
%     qvalue_ans = input(', wanna change it?[y/N] ', 's');
%     if strcmpi(qvalue_ans, 'y')
%         qvalue = input('? ');
%     end
% else
%     fprintf('qvalue: %s\n',qvalue)
% end
% 
% if imp_flag % Imputed variants --------------------------------------------
%     fprintf('Running cojo for %s\n', bedfile)
%     gcta_command = ".\gcta64.exe --bfile " + imp_dir + "\"+ ...
%         bedfile + " --cojo-file .\"+rawName+".pheno.ma"+...
%     " --cojo-slct --cojo-p "+ qvalue + " --cojo-wind 1000 --out "+...
%     rawName+"cojo"+bedRange+" --thread-num 25";
% %     gcta_command = ".\gcta64.exe --bfile " + imp_dir + "\"+ ...
% %         bedfile + " --cojo-file .\"+rawName+".pheno.ma"+...
% %     " --extract rs1.txt --cojo-cond rs2.txt --cojo-collinear 0.99 --cojo-p "+ qvalue + " --cojo-wind 1000 --out "+...
% %     rawName+"cojo"+bedRange+" --thread-num 25";
%     system(gcta_command);
%     delete(rawName+".pheno.ma")
%     if exist(rawName+"cojo"+bedRange+".jma.cojo", 'file')
%         cojo = readtable(rawName+"cojo"+bedRange+".jma.cojo", 'FileType', 'text');
%         save(rawName+"cojo"+bedRange+".mat", 'cojo')
%         delete(rawName+"cojo"+bedRange+".jma.cojo")
%         delete(rawName+"cojo"+bedRange+".cma.cojo")
%         delete(rawName+"cojo"+bedRange+".ldr.cojo")
%     end
%     delete(rawName+"cojo"+bedRange+".log")
%     fprintf('---------------------------------------------------------\n')
% 
% else % Genotyped variants -------------------------------------------------
%     gcta_command = ".\gcta64.exe --bfile ukb.qcn.cojo --cojo-file .\gcta.pheno.ma"+...
%         " --cojo-slct --cojo-p "+ qvalue + " --cojo-wind 1000 --out cojoout --thread-num 45";
%     system(gcta_command, '-echo');
%     delete('gcta.pheno.ma')
%     if exist('cojoout.cma.cojo', 'file')
%         delete("cojoout.cma.cojo")
%         delete("cojoout.ldr.cojo")
%     end
% end
% 
% end % END 
