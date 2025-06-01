function [geno, idx] = genoimpute(geno, opts)
% imputes and prunes genotype data.
% @03APR2023 Oveis Jamialahmadi, University of Gothenburg.

arguments
    geno
    opts.method {mustBeMember(opts.method, ["none", "random", "fixed", "bestguess"])} = "bestguess"
    opts.missingness (1,1) double = 0.15 % max missingness allowed
    opts.maf (1,1) double = 1 % remove MAF below this 
    opts.build_mask {mustBeMember(opts.build_mask, ["sum", "max", "comphet"])} = "max" % similar to REGENIE's build_mask argument
    opts.verbose (1,1) logical = true
end

if isstruct(geno)
    % bed field
    bim = geno;
    bim = rmfield(bim, "bed");
    geno = geno.bed;

    genostruct = true;
else
    genostruct = false;
end

geno = double(geno);
idx = isnan(geno) | (geno > 2) | (geno < 0); % missingness
geno(idx) = nan;
maf = mean(geno, 1, 'omitnan')./2;

% check missingness/MAF filters
rem_idx1 = sum(idx, 1)./size(idx, 1) > opts.missingness;
rem_idx2 = maf > opts.maf;
rem_idx =  rem_idx1 | rem_idx2;
if any(rem_idx)
    if opts.verbose
        if any(rem_idx1)
            fprintf('%d variants with a missingness > %.2g were removed.\n', ...
                sum(rem_idx), opts.missingness)
        end

        if any(rem_idx2)
            fprintf('%d variants with a MAF > %.2f were removed.\n', ...
                sum(rem_idx), opts.maf)
        end
    end
    geno(:, rem_idx) = [];
    
    if genostruct
        fis = string(fieldnames(bim.bim));
        for k = 1:numel(fis)
            bim.bim.(fis(k))(rem_idx) = [];
        end
    end
    
    idx = isnan(geno);
end

if opts.method == "none" % remove missigness
    idx = any(idx, 2);
    geno(idx, :) = [];
    
    geno = buildmask(geno, opts);
    if genostruct, bim.fam(idx, :) = []; end
    bim.bed = geno;
    geno = bim;
    clear bim
    return
end

p = size(geno, 2);
for i = 1:p
    if any(idx(:, i))
        if opts.method == "random"
            geno(idx(:, i), i) = binornd(2, maf(i), sum(idx(:, i)));
        elseif opts.method == "fixed"
            geno(idx(:, i), i) = maf(i).*2;
        elseif opts.method == "bestguess"
            geno(idx(:, i), i) = round(maf(i).*2);
        end
    end
end

idx = false(size(geno, 1), 1); % keep all genotypes
if any(isnan(geno), 'all')
    error('genoimpute:something went wrong!')
end

geno = buildmask(geno, opts);

% add other fileds to geno
if genostruct
    bim.bed = geno;
    geno = bim;
    clear bim
end


end % END

%% subfunctions ===========================================================
function geno = buildmask(geno, opts)

if isfield(opts, "build_mask")
    % use REGENIE build_mask argument
    if opts.build_mask == "max"
        geno = max(geno, [], 2);
    elseif opts.build_mask == "sum"
        geno = sum(geno, 2);
    elseif opts.build_mask == "comphet"
         geno = sum(geno, 2);
         geno(geno > 2) = 2;
    else
        error("genoimpute:unknown build_mask %s", opts.build_mask)
    end
end

end