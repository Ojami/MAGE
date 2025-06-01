function bgenix(bgen, opts)
% indexes bgen files uisng bgen library
% 
% Project: https://enkre.net/cgi-bin/code/bgen/dir?ci=trunk
% 
% get and build:
% 	wget http://code.enkre.net/bgen/tarball/release/bgen.tgz
% 	tar -xvzf bgen.tgz
% 	cd bgen.tgz
% 
% only keep ./apps/bgenix and compress the rest for a lighter version when
% transferring along with MAGE toolbox.
% 
% Oveis Jamialahmadi, University of Gothenburg, 1 March 2023.

arguments
    bgen {mustBeFile}
    opts.path {mustBeFolder} = fullfile(fileparts(which("bgenix.m")), "bgen.tgz", "build", "apps") % built apps from BGEN
    opts.verbose (1,1) logical = true
end

bgen = string(bgen);
if ~endsWith(opts.path, filesep), opts.path = opts.path + filesep; end
opts.path = makeWSLpath(opts.path);
if ispc, pre = "wsl "; else, pre = ""; end

if opts.verbose, disp("starting the indexing..."); end
for k = 1:numel(bgen)

    if opts.verbose, fprintf('%d (of %d)\n ', k, numel(bgen)); end

    [~, binfo] = fileattrib(bgen(k));
    bgen(k) = string(binfo.Name);

    cmd = opts.path + "bgenix -g " + makeWSLpath(bgen(k)) + " -index";
%     bgenix_command = bgenix_home + "bgenix -g " + bgenfile{ii} + " -incl-rsids rs370652263 > a.bgen";
%     bgenix_command = bgenix_home + "bgenix -help";
    if opts.verbose
        system(pre + cmd)
    else
        [~, ~] = system(pre + cmd);
    end
    
    if opts.verbose, fprintf("done\n"); end


end

end % END