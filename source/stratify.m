function idx = stratify(X, opts)
% A simple function to stratify a numeric vector based on quantiles/bins of
% data, which calls MATLAB built-in discretize function
% 
% Oveis Jamialahmadi, October 2024, University of Gothenburg.

arguments
    X {mustBeNumeric, mustBeVector} % vector of input values

    % only one of the followings should be used. 'q' overrides 'bins'.
    opts.q {mustBeVector, mustBeNumeric} % vector of quantiles (e.g. 0.5 for median)
    opts.bins {mustBeVector, mustBeNumeric} % bins for stratification

    % IncludedEdge argument of discretize function
    opts.side {mustBeMember(opts.side, ["left", "right"])} = "left"
end

assert(any(isfield(opts, ["q", "bins"])), "At least either of 'q' or 'bins' should be provided!")

if isfield(opts, "q")
    opts.vals = quantile(X, opts.q);
else
    opts.vals = opts.bins;
end

idx = discretize(X, [-inf, opts.vals, inf], IncludedEdge=opts.side);

end %END
