function [bgi_data, chr_len_no_zero, bgen_var_len] = bgireader(snps, opts)
% reads variant information from bgi (bgen index) files. Input must conform
% to input of gwasreader structure.

% Oveis Jamialahmadi, Sahlgrenska Academy, June 2021.
% 
% @05SEPT2024: 'table' option was added, if true returns the bgi info as a
% table.

arguments
    snps {mustBeA(snps, 'struct')}
    opts.bgenhome {mustBeFolder}
    opts.verbose (1,1) logical = true
    opts.parallel (1,1) logical = false
    opts.table (1,1) logical = false % return as table?
end

[~, geno] = getbulkgeno(strtrim(snps.snp), strtrim(snps.chr),...
    'home', opts.bgenhome, 'optionsOnly', true);

chrList = string(snps.chr);
chrom = unique(snps.chr);
if opts.verbose
    fprintf('Fetching input variants from bgi file...\n')
end

[bgi_data, chr_len_no_zero, bgen_var_len] = deal(cell(numel(chrom), 1));
tm = tic;
bgenhome = opts.bgenhome;
chrPos = geno.chrPos;
verbose = opts.verbose;
tmpsnps = cell(numel(chrom), 1);
for k = 1:numel(chrom)
    idx = chrList == chrom(k);
    tmpsnps{k, 1} = snps.snp(idx);
end

if opts.parallel
    astable = opts.table;
    parfor k = 1:numel(chrom)
        bgifile = fullfile(bgenhome, string(sprintf(chrPos, chrom(k)))) + ".bgi";
        if ~isfile(bgifile)
            if verbose, fprintf('%s does not exist!\n', bgifile), end
        end
        [bgi_data{k, 1}, chr_len_no_zero{k, 1}, bgen_var_len{k, 1}] = bgiloader(bgifile, tmpsnps{k}, verbose, astable);
    end
else
    for k = 1:numel(chrom)
        bgifile = fullfile(bgenhome, string(sprintf(chrPos, chrom(k)))) + ".bgi";
        if ~isfile(bgifile)
            if verbose, fprintf('%s does not exist!\n', bgifile), end
        end
        [bgi_data{k, 1}, chr_len_no_zero{k, 1}, bgen_var_len{k, 1}] = bgiloader(bgifile, tmpsnps{k}, verbose, opts.table);
    end
end

tm = toc(tm);
if opts.verbose
    fprintf('elapsed time: %.1f sec \n', tm)
end
end % END

%% subfunctions
function [bgi_data, chr_len_no_zero, bgen_var_len] = bgiloader(bgi_file, varInfo, verbose, astable)
conn = sqlite(bgi_file);
if ~isstring(varInfo)
    varInfo = string(varInfo);
end
% find position based queries: they must be in format of pos1-pos2
checkPos = arrayfun(@(x)split(x,'-'), varInfo, 'uni', false);
checkPosNumel = cellfun(@numel, checkPos); 
checkPosIdx1 = checkPosNumel == 2; % condition 1: numel must be 2
checkPosDouble = cellfun(@double, checkPos, 'uni', false);
checkPosIdx2 = cell2mat(cellfun(@(x)all(~isnan(x)), checkPosDouble, 'uni', false)); % condition 2: must be double
checkPosIdx = checkPosIdx1 & checkPosIdx2;

varInfo_rs = varInfo(~checkPosIdx);
varInfo_pos = checkPos(checkPosIdx);
START_str = "SELECT * FROM Variant WHERE ";
if ~isempty(varInfo_rs)
    varInfo_rs = "('" + join(varInfo_rs,"','") + "')";
    sqlquery_rs = "rsid IN "+varInfo_rs;
else
    sqlquery_rs = "";
end

if ~isempty(varInfo_pos)
    for k = 1:numel(varInfo_pos)
        if k > 1
            sqlquery_pos = sqlquery_pos + "(position BETWEEN " + varInfo_pos{k}(1) + " AND " + varInfo_pos{k}(end) + ")";
        else
            sqlquery_pos = "(position BETWEEN " + varInfo_pos{k}(1) + " AND " + varInfo_pos{k}(end) + ")";
        end
        if k < numel(varInfo_pos)
            sqlquery_pos = sqlquery_pos + " OR ";
        end
    end
else
    sqlquery_pos = "";
end
if isempty(varInfo_rs) || isempty(varInfo_pos)
    OR_str = "";
else
    OR_str = " OR ";
end
bgi_data = fetch(conn, START_str+sqlquery_rs+OR_str+sqlquery_pos);
close(conn)

if isempty(bgi_data) && verbose
    varInfo = cellstr(varInfo);
    error('ERROR: query cannot be found: \n%s\n', varInfo{:})
end

bgen_var_len = 0;
chr_len_no_zero = 0;

% Convert to struct
bgi_hdrs = {'chr', 'pos', 'snp', 'num_allele', 'a1', 'a2',...
    'file_start_position', 'size_in_bytes'};
try % < R2022a
    bgi_data = cell2struct(bgi_data', bgi_hdrs);
catch % >= R2022a
    bgi_data.Properties.VariableNames = bgi_hdrs;
    if ~astable
        bgi_data = table2struct(bgi_data, "ToScalar", false);
    end
end
end