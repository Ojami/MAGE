function [out1, out2] = interactionHelper(df, opts)
% interactionHelper performs interaction analysis for GLM (linear or
% logistic for now) using robust Sandwich HC3 estimator. This function
% calls directly interactionHelper.R function.
% Oveis Jamialahmadi, University of Gothenburg, May 24th 2025.

arguments
    df {mustBeA(df, "table")} % input table 
    opts.rfun {mustBeFile} = fullfile(fileparts(mfilename("fullpath")), "interactionHelper.R") % interactionHelper R code
    opts.dir {mustBeTextScalar} = pwd % dir to write/read the df to/from
    opts.formula_str {mustBeTextScalar} % formula to be passed to the R code
    opts.out {mustBeTextScalar} = "results" % output name prefix
    opts.bt (1,1) logical % is outcome a binary trait?
    opts.transform {mustBeMember(opts.transform, ["none", "log", "log1p", "log10", "irnt"])} = "irnt" % transformation function for continuous variables.
    opts.firth_engine {mustBeMember(opts.firth_engine, ["brglm2", "logistf"])} = "logistf"
    opts.nonlinear_E {mustBeTextScalar}
    opts.use_firth (1,1) logical = false
end

if opts.bt
    opts.bt = "bt";
else
    opts.bt = "qt";
end

% keep variables terms in formula
cols = split(string(opts.formula_str), ["~", "+", ":", "*"]);
df = rmmissing(df(:, cols));
n = height(df); % total sample size
[n_case, n_control] = deal(nan);
outcome = extractBefore(string(opts.formula_str), "~");
outcome = strip(string(outcome));

if opts.bt == "bt"
    outcome_vals = df.(outcome);
    n_case = nnz(outcome_vals);
    n_control = n - n_case;

elseif opts.transform ~= "none"
    df.(outcome) = feval(opts.transform, df.(outcome));
end

r = string;
r(1) = "setwd('" + replace(opts.dir, "\", "/") + "')";
r(numel(r) + 1) = "source('" + opts.rfun.replace("\", "/") + "')";

[df, tmp_file] = table2Rdf(df, write=true, ...
    dir=opts.dir, ...
    dt="parquet");

r = [r, df'];

ropts = opts;
fis = ["formula_str", "out", "bt", "firth_engine", "nonlinear_E", "use_firth"];
fis = setdiff(fieldnames(opts), fis);
ropts = rmfield(ropts, fis);
r(numel(r) + 1) = "interactionHelper(df, " + ...
    struct2rfun(ropts, replace=["out,prefix", "bt,outcome_type"]) + ")";

MATLAB2Rconnector(matlab.lang.makeValidName(opts.out) + ".r", ...
    code=r, delr=true);
delete(tmp_file)

% read the results
out1 = readtable(opts.out + ".summary.tsv", ...
    TextType="string", ...
    FileType="text",...
    VariableNamingRule="preserve");

out2 = readtable(opts.out + ".lh.tsv", ...
    TextType="string", ...
    FileType="text",...
    VariableNamingRule="preserve");

delete(opts.out + ".summary.tsv")
delete(opts.out + ".lh.tsv")

if opts.bt == "bt"
    out1 = renamevars(out1, ["Estimate", "Pr(>|z|)", ...
        "Std. Error", "Var1", "z value"], ...
        ["Beta", "P", "SE", "Term", "zValue"]);
    out1.zValue = [];
else
    out1 = renamevars(out1, ["Estimate", "Pr(>|t|)", ...
        "Std. Error", "t value", "Var1"], ...
        ["Beta", "P", "SE", "tValue", "Term"]);
    out1.tValue = [];
end

out1.N(:) = n;
out1.N_case(:) = n_case;
out1.N_control(:) = n_control;

end % END