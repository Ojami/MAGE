function pth = flashpca
% for copmiling: https://github.com/gabraham/flashpca 
% Note that spectra version should be < 1 (1 or higher doesn't work)
% Keep this function within flashpca directory where the executable file is
% Other files used to generate the executable file have been zipped (not
% need for running)
% 
% Oveis Jamialahmadi, 14NOV2024.

pth = fullfile(fileparts(which("flashpca.m")), "flashpca");

if ispc
    pth = makeWSLpath(pth);
end

end