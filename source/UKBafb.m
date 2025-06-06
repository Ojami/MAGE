function R = UKBafb(region)
% fetches info from The UK Biobank Allele Frequency Browser.
% 
% Oveis Jamialahmadi, September 2024, University of Gothenburg.

arguments
    region {mustBeTextScalar} % e.g.: "chr1-55039549-55039570"
end

getURL = "https://afb.ukbiobank.ac.uk/api";
headers = {'Content-Type' 'application/json'; 'Accept' 'application/json'};
options = weboptions('RequestMethod', 'post', 'HeaderFields', headers, ...
'Timeout', 20000);

data = struct;
data.query = ['query getRegion($regionStr: String!) {',newline,...
        '  region(regionStr: $regionStr) {', newline,...
        '    variants {', newline,...
        '      Chrom', newline,...
        '      Pos', newline,...
        '      Ref', newline,...
        '      Alt', newline,...
        '      rsID', newline,...
        '      maxImpact', newline,...
        '      alleleFreq', newline,...
        '      maxConsequence', newline,...
        '      geneSymbol}', newline,...
        '  }',newline,...
        '}'];
data.variables = struct("regionStr", region);

try
    R = webwrite(getURL, data, options);
    if isfield(R, "errors")
        disp(R.errors.message)
        return
    end
    R = R.data.region.variants;
    R = struct2table(R);
    R = convertvars(R, setdiff(colnames(R), ["Pos", "alleleFreq"]), @string);
catch ME
    R = ME.message;
end

end %END