function phetab = ldclump(assocfile, bedfile, opts)

% performs LD clumping using PLINK.
% @19MAY2023: realmin was changed to realmin*eps (1e-324)

arguments
    assocfile % file or a table
    bedfile {mustBeFile}
    opts.p1 (1,1) double {mustBeInRange(opts.p1, 0, 1)} = 5e-8
    opts.p2 (1,1) double {mustBeInRange(opts.p2, 0, 1)} = 0.01
    opts.r2 (1,1) double {mustBeInRange(opts.r2, 0, 1)} = 0.05
    opts.kb (1,1) double {mustBeGreaterThan(opts.kb, 0)} = 1000;
    opts.imp (1,1) logical = true % imputed genotypes
    opts.parallel (1,1) logical = false % only used for reading GWAS summary stat file
    opts.plinkdir {mustBeFolder} = fullfile(fileparts(which("ldclump.m")), "PLINK") % PLINK directory
    opts.verbose (1,1) logical = false % show PLINK output
    opts.save (1,1) logical = true % write clumped variants to a mat file
    opts.log10p (1,1) logical = false % P-values are in log10 scale. If true, the function converts them to P scale
end

if ~isfile(fullfile(opts.plinkdir, 'plink.exe'))
    error('couldn''t find plink.exe within %s!', opts.plinkdir)
end

if ~isstring(bedfile)
    bedfile = string(bedfile);
end
[wdir, bedfile, fileExt] = fileparts(bedfile);
bedfile = bedfile + fileExt;

%@30OCT2024: support for bgen files
if all(fileExt == ".bgen")
    opts.bgen = true;
else
    opts.bgen = false;
end

if ~opts.imp % not needed for imputed genotypes
    % Only include SNPs in the bed files
    bimfile = strrep(bedfile, '.bed', '');
    bimfile = bimfile + ".bim";
    getsnps = readtable(fullfile(wdir, bimfile), 'FileType', 'text');
    getsnps_var = getsnps.Properties.VariableNames;
    if numel(getsnps_var) > 1 % bim
        getsnps = getsnps.Var2;
    else % txt
        getsnps = getsnps.Var1;
    end   
end

% assocfile can be mat file (from BOLT-LMM/SAIGE/REGENIE) or PLINK2 --glm file
if istable(assocfile)
    phetab = assocfile;
    assocfile = getRandomName("ldclump", 5);
    toolFlag = 'MAT';
else
    [phetab, toolFlag] = readGWASfile(assocfile, 'full', true, ...
        'light', true, 'p', opts.p1, ...
        'parallel', opts.parallel);
end

if opts.log10p
    phetab.p = 10.^-phetab.p;
end
phetab(phetab.p > opts.p1, :) = [];

if isempty(phetab)
    fprintf('ldclump::no significant signal was found!\n')
    return
end
phetab.Properties.VariableNames(ismember(phetab.Properties.VariableNames, 'annotate')) = {'snp'};

% snpfix = phetab.snp;
% Change SNP names based on varid in BIM files
if ~opts.bgen
    [phetab, chrtag, bedfile] = snptag(phetab, fullfile(wdir, bedfile));
end

% Since p-value cannot go under realmin, modify 0 p-values based on
% the corresponding inverse Z-scores 
f_zero = phetab.p == 0; %@21MAY2025: strictly zero (user may have changed it before to realmin*eps)
if any(f_zero)
    Zinv = abs(phetab.se(f_zero)./phetab.beta(f_zero));
    phetab.p(f_zero) = 10.*eps*realmin.*(Zinv./max(Zinv));
end
%--------------------------------------------------------------------------
% Do LD-clumping using plink2cmd
% inCommands.plinkmethod = 1;

% p-value adjustment (mainly for BOLT-LMM output) -------------------------
% In BOLT-LMM, due to few digits used in fractional part, there are
% duplicates in p-values  which can be problematic for SNPs in LD; so,
% PLINK (probably) selects one at random as index. However, this may be not
% precise as their test statistics differ (as seen in COJO for instance). A
% better approach would be to use Z-scores to modify duplicate p-values to
% avoid this issue (only for nominally significant p-values). Still, it
% should be warned that BETA and SE of BOLT-LMM are for P_BOLT_LMM_INF and
% not P_BOLT_LMM (though the are highly correlated).
if contains(toolFlag, "BOLT"); phetab = adjdupP(phetab); end

% prepare the clump input file
phetab.Properties.VariableNames = upper(phetab.Properties.VariableNames); 
phetab.CHR = string(phetab.CHR);

% generate output names
[~, rawName] = fileparts(assocfile);

if opts.bgen
    chrtag = natsort(unique(phetab.CHR));
    rawName = rawName + "." + chrtag;
    genoPattern = strshare(bedfile, pattern=true);
    bedfile = compose(genoPattern, double(chrtag));
    bedfile = fullfile(wdir(1), regexprep(bedfile, '.bgen$', ''));
else

    rawName = rawName + "." + chrtag;
    bedfile = regexprep(bedfile, '.bed$', '');
end

genoline = strings(numel(bedfile), 1);
for k = 1:numel(bedfile)
    if opts.bgen
        genoline(k) = "--bgen " + bedfile(k) + ".bgen ref-first --sample " + bedfile(k) + ".sample";
    else
        genoline(k) = "--bfile " + bedfile(k) + "";
    end
end
    
%     inCommands.plinkcommand = "plink --threads 25 --bfile " + '"' + ...
%         fullfile(wdir,bedfile) + '"' + ...
%         " --r2 --ld-window-kb 1000 --ld-window 99999 --ld-snp-list snp.txt --out xx";

% swap ID and SNP cols to be consistent with BIM files
if ~opts.bgen
    SNP = phetab.SNP;
    phetab.SNP = phetab.ID;
    phetab.ID = [];
end

% begin the clumping
if opts.p2 < opts.p1
    opts = rmfield(opts, "p2");
end

clump = cell(numel(bedfile), 1);
for i = 1:numel(bedfile)
    
    if opts.verbose
        plinkTime = tic;
        fprintf('running clumping for file %d of %d...', i, numel(bedfile))
    end

    snpList = phetab.SNP(phetab.CHR == chrtag(i));
    if isempty(snpList) && opts.verbose
        fprintf('\b\b done!\n')
        fprintf('======================================\n')
        continue
    else
        writematrix(phetab.SNP(phetab.CHR == chrtag(i)), rawName(i) + ".txt",...
            'FileType', 'text', 'QuoteStrings', false)
    end
    
    writetable(phetab(phetab.CHR == chrtag(i), :),...
        rawName(i) + ".assoc", ...
        'FileType', 'Text', 'Delimiter', '\t')
    
    cmd = "plink --threads 45 " + ...
        genoline(i) + " --clump " + ...
        rawName(i) + ".assoc" + " --clump-p1 " + string(opts.p1) + ...
        " --clump-r2 " + opts.r2 + " --clump-kb " + opts.kb +...
        " --out " + rawName(i) + " --extract " + rawName(i) + ".txt";
    if isfield(opts, "p2")
         cmd = cmd + " --clump-p2 " + string(opts.p2);
    end

    if opts.bgen
        cmd = replace(cmd, textBoundary("start") + "plink", "plink2");
    end

    [~, ~] = system(opts.plinkdir + string(filesep) + cmd);

    if opts.verbose
        fprintf('\b\b done!\n');
        toc(plinkTime)
        fprintf('======================================\n')
    end
    
    warning off
    delete(rawName(i) + ".txt")
    delete(rawName(i)+".log")
    delete(rawName(i) + ".assoc")
    warning on
    
    % Read clumped file
    if ~exist(rawName(i)+".clumped", 'file')
        % No significant clumping results.
        continue
    end

    clumptab = readtable(rawName(i)+".clumped",...
        'FileType', 'text', 'ReadVariableNames',...
        true, 'format', '%f%f%s%f%f%f%f%f%f%f%f%s');
    
    warning off
    delete(rawName(i)+".clumped")
    warning on
    clump{i} = string(clumptab.SNP);
end

clump(cellfun(@isempty, clump)) = [];
clump = vertcat(clump{:});

if isempty(clump)
    if opts.verbose, fprintf('No independent SNP found!\n'); end
    return
end

if ~opts.bgen
    phetab.ID = phetab.SNP;
    phetab.SNP = SNP;
    phetab(~ismember(phetab.ID, clump), :) = [];
end

if opts.save
    [~, rawName] = fileparts(assocfile);
    assocfile = rawName + ".clump.mat";
    save(assocfile, 'phetab')
end

if opts.verbose
    fprintf('%d lead SNPs were written to %s\n', numel(clump), assocfile)
end
end

%% subfunctions ===========================================================
function [subgss, chr, bedfiles] = snptag(gss, bimfiles)

bimfiles = regexprep(bimfiles, '.bed$', '.bim');
if ~isstring(bimfiles)
    bimfiles = string(bimfiles);
end

% get chromosomes in bim files
chr = strings(numel(bimfiles), 1);
for i = 1:numel(bimfiles)
    fid = fopen(bimfiles(i), 'r');
    bimLine = string(split(fgetl(fid)));
    chr(i) = bimLine(1);
    fclose(fid);
end

% only keep summary stats with these chromosomes
gss(~ismember(string(gss.chr), chr), :) = [];

% only keep required bim files
rmidx = ~ismember(chr, string(unique(gss.chr)));
bimfiles(rmidx) = []; chr(rmidx) = [];
bedfiles = regexprep(bimfiles, '.bim$', '.bed');

if ~any(colnames(gss) == "allele0")
    % @29MAY2025: the logic of treating a2 as allele1 (ergo effect allele)
    % is not consistent with most tools. Therefore a2 is treated as
    % alternate allele from now on.
    try
        gss = renamevars(gss, ["a2", "a1"], ["allele0", "allele1"]);
        fprintf("\tldclump: a1(allele1) is effect and a2(allele0) alternate allele.\n")
    catch
         gss = renamevars(gss, "allele2", "allele0");
         fprintf("\tldclump: allele1 is effect and allele2(allele0) is alternate allele.\n")
    end
    fprintf("\t\t\tif these do not correspond to your summary stats, rerun with appropriate col names\n")
end
gss.id = join([gss.chr, gss.pos, gss.allele0, gss.allele1], ":");

chrstr = string(gss.chr);
subgss = cell(numel(bimfiles), 1);
for ii = 1:numel(bimfiles)

    cidx = chrstr == chr(ii);
    subgss{ii} = gss(cidx, :);
    gss(cidx, :) = []; chrstr(cidx) = [];

    % read bim file (only those matching the variants in the input table)
    ids_extract = [subgss{ii}.chr + ":" + ...
        subgss{ii}.pos + ":" + ...
        subgss{ii}.allele0 + ":" + ...
        subgss{ii}.allele1; ...
        subgss{ii}.chr + ":" + ...
        subgss{ii}.pos + ":" + ...
        subgss{ii}.allele1 + ":" + ...
        subgss{ii}.allele0;...
        subgss{ii}.snp];
    bim = bimreader(bimfiles(ii), struct=false, ...
        cols=["chr", "snp", "cmg", "pos", "a1", "a2"], ...
        extract=ids_extract);
    clear ids_extract

    % if they match, return
    if all(ismember(subgss{ii}.snp, bim.snp))
        subgss{ii}.id = subgss{ii}.snp;
        continue
    end
    
    %@16MAY2025: this approach is deprecated, instead throw an error
    % check if there exists duplicate (multi-allelic) variants. If so,
    % modify original bim files (backed up already).
    dups = duplicates(bim.snp);
    if ~isempty(dups)

        error("duplicated snps in BIM files are not allowed!")

        % % modify original bim file
        % didx = ismember(bim.snp, dups);
        % bim.snp(didx) = join([bim.chr(didx), bim.pos(didx), ...
        %     bim.a1(didx), bim.a2(didx)], ":");
        % writetable(bim, bimfiles(ii), "WriteVariableNames", false, ...
        %     "Delimiter", char(9), "FileType", "text")
    end

    bim.id = join([bim.chr, bim.pos, bim.a1, bim.a2], ":");

    % check non-matching varinats
    midx = ~ismember(subgss{ii}.id, bim.id);
    if any(midx)
        % try flipping alleles and check again
        subgss{ii}.idf = join([subgss{ii}.chr, subgss{ii}.pos, subgss{ii}.allele1, subgss{ii}.allele0], ":");
        mfidx = ~ismember(subgss{ii}.idf, bim.id) & midx;
        if any(mfidx) % cannot be found within bim files, should be removed from gss
            subgss{ii}(mfidx, :) = [];
        end
        
        % replace variant ids with flipped ones for these variants
        fidx = ismember(subgss{ii}.idf, bim.id);
        subgss{ii}.id(fidx) = subgss{ii}.idf(fidx);
        subgss{ii}.idf = [];
    end
    
    [idx1, idx2] = ismember(subgss{ii}.id, bim.id);
    assert(all(idx1), "something went wrong!")
    subgss{ii}.id = bim.snp(idx2);

    % [idx1, idx2] = ismember(bim.id, subgss{ii}.id);
    % subgss{ii}.id(idx2(idx1)) = bim.snp(idx1);

end
subgss = vertcat(subgss{:});

end
%%
function phetab = adjdupP(phetab)
% modifies duplicate p-values based on BETA/SE (Z-scores).
% INPUTS:
%   phetab: a table containing summary stat with following columns
%   required: chr, p, beta, se

chrLoop = unique(phetab.chr);
adjPtop = -1.*ones(size(phetab, 1), 1);
for i = 1:numel(chrLoop) % loop over chr
    chrIdx = phetab.chr == chrLoop(i);
    if sum(chrIdx) ~= numel(unique(phetab.p(chrIdx))) % any duplicates in chr?
        % create a temp table 
        ptab = ismember(phetab.Properties.VariableNames, {'p', 'beta', 'se'});
        ptab = phetab(chrIdx, ptab);
        adjP = -1.*ones(size(ptab, 1), 1);  
        Pcheck = ptab.p <= 0.05; % check all nominally significant signals

        if any(Pcheck) % any significant signal for this CHR
            if sum(Pcheck) ~= numel(unique(ptab.p(Pcheck))) % any duplicate in CHR?
                Pcheck = find(Pcheck);

                [Pu, ~, Pi] = unique(ptab.p(Pcheck));
                Pcount = histcounts(Pi,[1:numel(Pu),inf]);
                Pu = Pu(Pcount > 1);  
                if ~isempty(Pu) % found duplicates
                    f_dupCheck = [];
                    for j = 1:numel(Pu)
                        f_dup = Pcheck(ptab.p(Pcheck) == Pu(j));
                        Z = abs(ptab.beta(f_dup)./ptab.se(f_dup)); % Z-score
                        [~, Z] = sort(Z, 'descend');
                        adjP(f_dup) = ptab.p(f_dup).*(1 + Z./1e5);
                        if numel(unique(adjP(f_dup))) ~= numel(adjP(f_dup))
                            erro('something went wrong with unique!')
                        end
                        f_dupCheck = [f_dupCheck; f_dup];
                    end
                end
            end
        end
        
        % check if order of adjusted p-values is preserved within max margins
        porderFlag = checkOrder(phetab.p(chrIdx), adjP, f_dupCheck);
        if ~porderFlag
            error('p-value adjustment changed order of p-values in other clumps.')
        end
        
        lastCheckP = adjP(Pcheck);
        lastCheckP(lastCheckP == -1) = [];
        if numel(unique(lastCheckP)) ~= numel(lastCheckP)
            error('duplicate found after last check-up!')
        end
    
        adjPtop(chrIdx) = adjP;
    end
   
end

if any(adjPtop ~= -1)
    fprintf('%d p-value were adjusted.\n', sum(adjPtop ~= -1))
    phetab.p(adjPtop ~= -1) = adjPtop(adjPtop ~= -1);
else
    fprintf('no p-value adjustment is required.\n')
end
end

%%
function orderFlag = checkOrder(ptab, adjP, expIdx)
Padjusted = ptab;
Padjusted(adjP ~= -1) = adjP(adjP ~= -1);
[~, Poriginal] = sort(ptab);
[~, PadjIdx] = sort(Padjusted);
discrep = Poriginal(Poriginal ~= PadjIdx);
if all(ismember(discrep, expIdx))
    orderFlag = true;
else
    orderFlag = false;
end
end

