function runVEP(infiles, opts)
% annotates input (e.g. vcf) files using ensembl-vep in 'pth' path. Only
% supports homo_sapiens now; this can be modified by adding appropriate
% flags/options as described in:
% http://www.ensembl.org/info/docs/tools/vep/script/vep_options.html
% 
% Oveis Jamialahmadi, Sahlgrenska Academy, University of Gothenburg, March 2022.

arguments
    infiles {mustBeFile, mustBeVector}
    opts.pth = "" % leave empty if already added to path "/home/ensembl-vep/"
    opts.dir = "$HOME/.vep/" % "E:\VEP\ensembl-vep\cache\"
    opts.fasta {mustBeTextScalar} % path to fasta file. when cache/offline is active it should automatically detect it. But if returned error, uncomment following lines for fasta and 
    % make sure Homo_sapiens.GRCh[37|38].dna.primary_assembly.fa file is there
    opts.species = "homo_sapiens"
    opts.assembly {mustBeMember(opts.assembly, ["GRCh37", "GRCh38"])} = "GRCh38" % for --species homo_sapiens
    opts.cache (1,1) logical = true % cache/fasta files must be present in $HOME/.vep
    opts.offline (1,1) logical = true
    opts.force_overwrite (1,1) logical = true
    opts.symbol (1,1) logical = true
    opts.check_existing (1,1) logical = true
    opts.fork (1,1) double = 30 % number of parallel jobs 
    opts.no_stats (1,1) logical = true
    opts.check_ref (1,1) logical = true % needs --fasta (by default VEP finds it in $HOME/.vep)
    opts.dont_skip (1,1) logical = false
    opts.buffer_size (1,1) double = 1.5e5 % 5e3 % variant chunk number. VEP default is 5e3
    opts.tab (1,1) logical = true 
    opts.nearest {mustBeTextScalar, mustBeMember(opts.nearest, ["none", "transcript", "gene", "symbol"])} = "none"
    opts.biotype (1,1) logical = false
    opts.domains (1,1) logical = false
    opts.merged (1,1) logical = false % it true, uses _merged cache files (Ensembl + RefSeq)
    opts.refseq (1,1) logical = false % if ture, uses _refseq cache files inconsistent with 'merged'
    opts.most_severe (1,1) logical = false
    opts.fields {mustBeText, mustBeVector} = ["Uploaded_variation", ...
        "Location", "Gene", "Consequence", "Existing_variation", "SYMBOL", ...
        "IMPACT", "Amino_acids", "CHECK_REF", "DOMAINS"] % BIOTYPE is for transcript and not protein 

    % non-native flags
    opts.parse (1,1) logical = true % to parse VEP output 
    opts.parallel (1,1) logical = true % to process tall tables
    opts.workers (1,1) double = 40
    opts.merge (1,1) logical = true % merge all single files (only if numel infiles > 1), and write into a mat and txt files
end

% check inputs
if contains(opts.pth, "\"); opts.pth = makeWSLpath(opts.pth); end
if ~isstring(infiles); infiles = string(infiles); end
if ispc
    opts.dir = makeWSLpath(opts.dir);                                    
end

% if ~isfield(opts, 'fasta') && opts.check_ref
%     opts.fasta = opts.dir;
% end

chcmd = ""; % "sudo chmod -R a+rwx " + opts.dir;
% if ispc % or 777
%     % system("wsl sudo chmod -R a+rwx " + opts.dir);
% else
%     % system("sudo chmod -R a+rwx " + opts.dir);
% end

opts.dir = '"' + opts.dir + '"';
% if isfield(opts, 'fasta')
%     opts.fasta = '"' + opts.fasta + '"';
% end

if opts.most_severe % remove inconsistent flags with most_severe
    most_severe_to_flip = ["appris", "biotype", "canonical", "ccds", ...
        "coding_only", "domains", "flag_pick", "flag_pick_allele", ...
        "no_intergenic", "numbers", "pick", "pick_allele", "polyphen", ...
        "protein", "sift", "summary", "symbol", "tsl", "uniprot", ...
        "xref_refseq"];
    opts = nullifyVEPflags(opts, most_severe_to_flip);
end
opts.fields = join(opts.fields, ",");
if ~opts.check_ref
    opts.fields(opts.fields == "CHECK_REF") = [];
end

% check parallel pool to parse VEP resuls
if ~opts.parse; opts.parallel = false; end
if opts.parallel && isempty(gcp('nocreate'))
    parpool('Processes', opts.workers);
end

% construct the vep cmd pattern
rmfis = ["parse", "parallel", "merge", "workers"];
vep = parseVEPcmd(opts, rmfis);

for k = 1:numel(infiles)
    t1 = tic;
    fprintf('%d-annotating %s (of %d) ', k, infiles(k), numel(infiles))
    outfile = regexprep(infiles(k), '.vcf$', '.txt');
    if ispc, dos2unix(makeWSLpath(infiles(k))); end
    vepcmd = compose(vep, makeWSLpath(infiles(k)), makeWSLpath(outfile));
    if ~isfile(outfile)
        runbash([chcmd; vepcmd], "vep", "parallel", false, "wait", true, "verbose", false)
    end

    if opts.parse
        fprintf('\n\tparsing...')
        if ~isfile(regexprep(outfile, '.txt$', '.mat'))
            vepParser(outfile, opts)
        end
        fprintf('\b\b')
    end

    telapsed = seconds(toc(t1));
    telapsed.Format = "hh:mm:ss";
    fprintf("(done in %s)\n", telapsed)
end

if opts.merge && numel(infiles) > 1 % merge all
    tic
    output = "vep.merged";
    matfiles = regexprep(infiles, '.vcf$', '.mat');
    matfiles(~isfile(matfiles)) = [];
    if isempty(matfiles)
        disp("no file has remained to merge!")
        return
    end
    matfiles = natsort(matfiles);
    tab = fileDatastore(matfiles,'ReadFcn', @(x) load(x).tab);
    tab = gather(tall(tab));
    tab = vertcat(tab{:});
    tab.Properties.VariableNames = erase(tab.Properties.VariableNames, '#');
    fprintf('writing to %s ...', output + ".txt")
    writetable(tab, output + ".txt", 'Delimiter', "\t", ...
        "WriteVariableNames", true)
    fprintf('\b\b Done.\n')
    
    fprintf('writing to %s ...', output + ".mat")
    catcols = ["Gene", "Consequence", "SYMBOL", "IMPACT"];
    if opts.domains
        catcols = [catcols, "DOMAINS"];
    end
    tab = convertvars(tab, catcols, "categorical");
    tab = table2struct(tab, 'ToScalar', true);
    save(output + ".mat", "-struct", 'tab')
    fprintf('\b\b Done.\n')
    toc
    
end

end % END

%% subfunctions ===========================================================
function opts = nullifyVEPflags(opts, nullFlags)
fnames = string(intersect(fieldnames(opts), nullFlags));
for i = 1:numel(fnames)
    if islogical(opts.(fnames(i)))
        opts.(fnames(i)) = false;
    else
        opts.(fnames(i)) = "none";
    end
end

end

%% ------------------------------------------------------------------------
function cmd = parseVEPcmd(opts, rmfis)
pth = opts.pth;
opts = rmfield(opts, union("pth", rmfis));
if ~endsWith(pth, "/") && pth ~= ""; pth = pth + "/"; end
cmd = pth + "vep -i %s -o %s";

fnames = string(fieldnames(opts));
for i = 1:numel(fnames)
    if islogical(opts.(fnames(i))) 
        if opts.(fnames(i)); cmd = cmd + " --" + fnames(i); end
    elseif string(opts.(fnames(i))) ~= "none"
        if contains(string(opts.(fnames(i))), ",")
            cmd = cmd + " --" + fnames(i) + ' "' + opts.(fnames(i)) + '"';
        else
            cmd = cmd + " --" + fnames(i) + " " + opts.(fnames(i));
        end
    end
end

end

%% ------------------------------------------------------------------------
function vepParser(outfile, opts)

% HIGH/MODERATE impacts from http://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences
filters = ["transcript_ablation", "splice_acceptor_variant",...
    "splice_donor_variant", "stop_gained", "frameshift_variant", ...
    "stop_lost", "start_lost", "transcript_amplification", ...
    "inframe_insertion", "inframe_deletion", "missense_variant"];
impacts = ["HIGH", "MODERATE"];
matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.NullProgressReporter);
tab = tabularTextDatastore(outfile, 'TextType', 'string', ...
    'CommentStyle', '##', 'NumHeaderLines', 0, ...
    'VariableNamingRule', 'preserve');

% remove unnecessary columns 
if opts.domains
    unneeded_cols = ["BIOTYPE", "CHECK_REF"];
else
    unneeded_cols = ["BIOTYPE", "CHECK_REF", "DOMAINS"];
end
tab.SelectedVariableNames = setdiff(tab.SelectedVariableNames, ...
    unneeded_cols, 'stable');
tab = tall(tab);

rmIdx = ~contains(tab.IMPACT, impacts);
tab(rmIdx, :) = [];
tab = gather(tab);
matlab.bigdata.internal.executor.ProgressReporter.override(matlab.bigdata.internal.executor.CommandWindowProgressReporter);
tab = unique(tab, 'rows');
tab.Properties.VariableNames = erase(tab.Properties.VariableNames, '#');

if isempty(tab), return, end

% keep only most severe consequence 
cons = arrayfun(@(x)unique(split(x, ",")), unique(tab.Consequence), 'uni', false);
cons = unique(vertcat(cons{:}));
cons = setdiff(cons, filters); % remove these consequences
tab.Consequence = erase(tab.Consequence, cons);
tab.Consequence = regexprep(tab.Consequence, ",+$|^,+", "");
tab(ismember(tab.Consequence, ["", ","]), :) = []; % e.g. protein_altering_variant has a moderate impact but not in desired filters
tab = unique(tab, 'rows');

% @22JUNE2023: we no longer merge overlapping variants since this may cause
% wrong (although only for few variants) functional (e.g. LoF) assignments
% at later steps. As example "1:151067336:T:C" variant overlaps CDC42SE1
% and MLLT11, one is splice_acceptor_variant and the other is
% missense_variant. So, it creates a dilemma for later tools to dissect
% which is which (are both LoF if the tool decides to keep the most severe
% one? no!). Hence, we leave them as they are (i.e. those variants are
% repeated in the final file, once for each gene). However, we merge
% multiple consequences (due to VEP annotations from different trnascripts)
% for each variant.
tab.tmpID = tab.Uploaded_variation + ":" + tab.Gene;
[dups, ~, counts, tidx] = duplicates(tab.tmpID);
dupidx = ismember(tab.tmpID, dups);
tab_dup = groupsummary(tab(dupidx, :), 'tmpID', @(x)join(unique(x, 'stable'), ","));

% % merge overlapping varints (e.g. variants overlapping multiple genes).
% % These variants will be seperated in saigeWrapper/getWES for
% % gene/region-tests
% [dups, ~, counts, tidx] = duplicates(tab.Uploaded_variation);
% dupidx = ismember(tab.Uploaded_variation, dups);
% tab_dup = groupsummary(tab(dupidx, :), 'Uploaded_variation', @(x)join(unique(x, 'stable'), ";"));

tab_dup(:, ["tmpID", "GroupCount"]) = [];
tab.tmpID = [];
tab_dup.Properties.VariableNames = erase(tab_dup.Properties.VariableNames, "fun1_");
assert(all(colnames(tab) == colnames(tab_dup)))
tab = tab(tidx, :);
tab(counts > 1, :) = tab_dup;
save(regexprep(outfile, '.txt$', '.mat'), 'tab')
end
