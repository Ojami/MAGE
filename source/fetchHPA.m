function [tab, failed_proteins] = fetchHPA(plist, opts)
% Reads data from Human Protein Atlas website for a list of
% query proteins. Returns a table of HPA data for query proteins found in
% the HPA database.

% Oveis Jamialahmadi, University of Gothenburg. Jan 2024

arguments
    plist {mustBeText, mustBeVector} % list of input proteins. Can be a mixture of ENSG ids and official gene names.
    opts.gtf {mustBeFile} = fullfile(fileparts(which("fetchHPA.m")), "GTF", "Homo_sapiens.GRCh38.107.gtf.gene.mat") % GTF file if there are gene name in the query list
    opts.verbose (1, 1) logical = true
    opts.pause (1,1) double = 0 % pausing time in second between each fetch to avoid server overload (default is 0: no pause)
end

plist = string(plist);
plist = rmmissing(plist);

tab = table; 
failed_proteins = string;
if isempty(plist), return; end

% check if the list contains gene names instead of gene ID
idx = ~contains(plist, "ENSG" + digitsPattern);

if any(idx) % map gene names to gene IDs
    name_list = plist(idx);
    gtf = load(opts.gtf, "gene_name", "gene_id");

    [f1, f2] = ismember(name_list, gtf.gene_name);
    id_list = strings(numel(name_list), 1);
    id_list(:) = missing;
    id_list(f1) = gtf.gene_id(f2(f1));
    
    plist(idx) = id_list;
    plist = rmmissing(plist);
end

plist = unique(plist);

% initialize web options
headers = {'Content-Type' 'application/xml'};
woptions = weboptions('HeaderFields', headers, 'Timeout', 5000);

tmp = cell(numel(plist), 1);

t1 = tic;
if opts.verbose
    fprintf("Fetching data for %d proteins ", numel(plist))
    progressBar = progressGen(numel(plist)-1);
end

failed_fetch = zeros(numel(plist), 1);
for k = 1:numel(plist)
    try
        tmp{k} = readtable("http://www.proteinatlas.org/" + plist(k) + ".tsv", ...
            "TextType", "string", "FileType", "text", ...
            "VariableNamingRule", "preserve", "WebOptions", woptions);
    catch ME
        if startsWith(ME.message, "File not found at")
            failed_fetch(k) = 1;
        else
            failed_fetch(k) = 2;
        end

        continue
    end

    if opts.pause > 0, pause(opts.pause); end
    
    if opts.verbose, progressGen(progressBar, k); end
end

tab = vertcat(tmp{:});

% final summary
if opts.verbose
    not_found_pt = failed_fetch == 1;
    if any(not_found_pt)
        fprintf("\n%d proteins were not found in HPA: %s\n", ...
            sum(not_found_pt), ...
            join(plist(not_found_pt), ","))
    end

    er_pt = failed_fetch == 2;
    if any(er_pt)
        fprintf("\nfailed fetching from server due to other reasons for %n proteins: %s\n", ...
            sum(er_pt), join(plist(er_pt), ","))
    end

    if all(~failed_fetch) % nothing failed
        fprintf("\n") % for final time elapsed
    end
end

failed_proteins = plist(failed_fetch ~= 0);
t1 = seconds(toc(t1));
fprintf('Finished in %.2f sec\n', seconds(t1))

end % END