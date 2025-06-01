function pth = qctool(opts)

% returns absolute path to qctool binary file
arguments
    opts.path {mustBeFolder} = fileparts(which("qctool.m"))
    opts.unix (1,1) logical = true
end

fis = getfilenames(opts.path).x_;
fis = fis(startsWith(fis.lower, "qctool_v"));
fis = fis(1); % multiple versions?
pth = fullfile(opts.path, fis);

if opts.unix
    pth = makeWSLpath(pth);
end

end