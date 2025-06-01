function res = fetchReplicationGWAS(snps, ssfiles, opts)
% this function fetches summary stats from external GWAS summary stat
% files. These files can be from different resources, e.g. GWAS Catalogue,
% etc.
% 
% Oveis Jamialahmadi, University of Gothenburg, May 2023.

arguments
    snps {mustBeA(snps, "struct")}
    ssfiles {mustBeFile} % summary stat files
    opts.bgenhome {mustBeFolder} = "D:\Imputed" % to add 'pos' to snps
    opts.liftover (1,1) logical = true % map GRCh37 to 38 (GWAS data can be in GRCh38). 
end

if ~isfield(snps, "pos")
    pos = bgireader(snps, bgenhome=opts.bgenhome);
    pos = struct2table(vertcat(pos{:}));
    [~, idx] = ismember(snps.snp, pos.snp);
    snps.pos = pos.pos(idx);
end

fis = string(fieldnames(snps));
for k = 1:numel(fis)
    if isrow(snps.(fis(k)))
        snps.(fis(k)) = snps.(fis(k))';
    end
end

if opts.liftover
    li = liftover(struct2table(snps));
    [~, idx] = ismember(snps.snp, li.snp);
    snps.pos38 = li.pos(idx);
end
snps.patt = "(?<!-)\b" + snps.snp + "(?!-)\b";

[wd, name, ext] = fileparts(ssfiles);
wd(wd == "") = pwd;
ssfiles = name + ext;

out = cell(numel(snps.snp), 1);
for j = 1:numel(snps.snp)
    fprintf("variant %d of %d\n", j, numel(snps.snp))
    s.snp = snps.snp(j);
    s.chr = snps.chr(j);
    s.pos = snps.pos(j);
    if opts.liftover
        s.pos38 = snps.pos38(j);
    end
    s.patt = snps.patt(j);

    [gc, hdrs] = deal(cell(numel(ssfiles), 1));
    for k = 1:numel(ssfiles)

        fprintf("\t%d/%d-reading %s\n", k, numel(ssfiles), ssfiles(k))
        smr = bfilereader(fullfile(wd(k), ssfiles(k)), "summary", "only", "header", true, "return", "raw");
        smr = smr(1:2);
        smr = split(smr);

        hdr = smr(1, :);        
        idx = find(contains(hdr.lower, ["snp", "variant", "rsname", "rsid"]));
        idx_chr = find(startsWith(hdr.lower, "chr"));
        idx_pos = find(ismember(hdr.lower, ["bp", "position", "pos", "base_pair_location"]) | contains(hdr.lower, "genpos"), 1);
        
        if isempty(idx)
            tab = "";
        else
            tab = bfilereader(fullfile(wd(k), ssfiles(k)), "pattern", s.patt, ...
                "patternCol", idx, "return", "raw", "parallel", true);
            tab = split(string(tab'));
        end
    
        if all(tab == "")
            if isempty(idx_chr) || isempty(idx_pos), continue; end
            
            smr = smr(2, idx_chr);
            if smr.startsWith("chrom")
                s_chr = "chrom" + s.chr;
            elseif smr.startsWith("chr")
                s_chr = "chr" + s.chr;
            end

            tab1 = bfilereader(fullfile(wd(k), ssfiles(k)), ...
                "pattern", [string(s_chr), string(s.pos)], ...
                "patternCol", [idx_chr, idx_pos], "return", "raw", ...
                "parallel", true, "multiCol", true);
            tab1 = split(string(tab1'));
            
            % @25JULY2023: since solely based on summary stat name it's not
            % clear if positional are based on which reference genome, we
            % fetch both based on GRCh37 and 38, then later, the user may
            % further check which one is correct and keep only those
            % summary stats.
            % if all(tab == "") && opts.liftover 
            if opts.liftover % check also with GRCh38
                tab2 = bfilereader(fullfile(wd(k), ssfiles(k)), ...
                    "pattern", [string(s_chr), string(s.pos38)], ...
                    "patternCol", [idx_chr, idx_pos], "return", "raw", ...
                    "parallel", true, "multiCol", true);
                tab2 = split(string(tab2'));
                if all(tab2 == "")
                    tab = tab1;
                elseif all(tab1 == "")
                    tab = tab2;
                else
                    tab = [tab1; tab2];
                end
            else
                tab = tab1;
            end
        end
        
        if all(tab == ""), continue; end
        % tab = reshape(tab, numel(hdr), [])';
        if iscolumn(tab), tab = tab'; end
        hdrs{k} = createGWASheader(hdr);
        gc{k} = array2table(reshape(tab, numel(hdrs{k}), [])', VariableNames=hdrs{k});
        gc{k}.Pheno(:) = ssfiles(k).extractBefore(".");
        hdrs{k} = colnames(gc{k});
        if any(hdrs{k} == "OR") && ~any(hdrs{k} == "BETA")
            gc{k}.BETA = string(log(double(gc{k}.OR)));
            gc{k}.OR = [];
            hdrs{k} = colnames(gc{k});
        end
    end

    idx = cellfun(@isempty, gc);
    gc(idx) = []; hdrs(idx) = [];
    if isempty(hdrs)
        out{j} = [];
        continue
    end

    hdr = hdrs{1};
    for k = 2:numel(hdrs), hdr = union(hdr, hdrs{k}); end
    
    res = cell(numel(gc), 1);
    for k = 1:numel(gc)
        res_tmp = strings(height(gc{k}), numel(hdr));
        [f1, f2] = ismember(hdr, colnames(gc{k}));
        res_tmp(:, f1) = gc{k}{:, f2(f1)};
        res{k} = res_tmp;
    end
    res = array2table(vertcat(res{:}), VariableNames=hdr);
    res.SNP(:) = s.snp;
    out{j} = res;
end
out(cellfun(@isempty, out)) = [];

hdr = colnames(out{1});
for k = 2:numel(out), hdr = union(hdr, colnames(out{k})); end
res = cell(numel(out), 1);
for k = 1:numel(out)
    [f1, f2] = ismember(hdr, colnames(out{k}));
    tres = strings(height(out{k}), numel(hdr));
    tres(:, f1) = out{k}{:, f2(f1)};
    res{k} = tres;
end
res = array2table(vertcat(res{:}), VariableNames=hdr);
res = convertvars(res, ["P", "BETA", "SE"], @double);

end % END