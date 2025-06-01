function parseEpirr
% reads epigenetic meta-data from EiRR API for RoadMap project's cell names
% and adds shorter names to RoadMap_ENTEX_Finuacane2018_supp7.mat table for
% visualization purposes

tab = load("RoadMap_ENTEX_Finuacane2018_supp7.mat").tab;

list.cell = unique(tab.(1)); % unique cell-types
list.length = strlength(list.cell);
list = struct2table(list);
list(list.length <= 30, :) = [];
% list.clean = extractAfter(list.cell, "_");
list.clean = list.cell;
list.clean = replace(list.clean, "_", " ");
list.clean = strtrim(regexprep(list.clean, "skin\w+", ""));

% wopts = weboptions(Timeout=1e5);
% epir = webread("https://www.ebi.ac.uk/epirr/api/v1/epigenomes", wopts);
epir = load("epir.mat").epir;
epir = struct2table(epir);
epir = convertvars(epir, 1:width(epir)-1, @string);

for k = 1:height(list)
    term = list.clean(k).lower;
    idx = contains(epir.description.lower, term);
    if ~any(idx)
        idx = contains(epir.description.lower, strtrim(regexprep(term, "skin\w+", ""))); 
    end
end

end