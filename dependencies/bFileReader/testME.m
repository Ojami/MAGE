function testME
% clc
%% Add the folder containing bFileReaderDep.py to the Python path
py.sys.path().append('D:\MATLAB.dependencies\MAGE\bFileReader');

%% Import and reload the module
mod = py.importlib.import_module('bFileReaderDep');
mod = py.importlib.reload(mod);

%% Set common parameters
% Specify the file (uncompressed for testing; if using compressed, the code auto-detects)
filename = 'Cirrhosis_MVP.txt';  
% Use auto-detection for delimiter by passing py.None (or you can explicitly specify, e.g. ',' for CSV)
delim = py.None;  
skip = int32(0);
bench = 'sp';
blocksize = int32(5000);  % Use 1 MB blocksize to force partitioning

%% Test 1: Filter rows by SNP_ID using compareToCols
% This command filters rows where the value in column 0 (SNP_ID) exactly matches "^rs531691745$"
result_py = mod.bFileReaderDep.compareToCols(...
    filename, ...
    py.list({'^rs738409$'}), ...       % regex pattern list for SNP_ID
    delim, ...
    py.list({int32(0)}), ...              % column index 0 (SNP_ID)
    skip, ...
    bench, ...
    py.None, ...                        % retrunCols = None means return the full line
    blocksize);
result1 = char(result_py);
disp('Result from compareToCols (filtering SNP_ID):');
disp(result1);

% %% Test 2: Filter rows by chrom and pos using filterCol
% % Here we want rows where chrom == 1 and pos == 284844.
% % Given the header (0-indexed): SNP_ID=0, chrom=1, pos=2, we filter on columns 1 and 2.
% result_py = mod.bFileReaderDep.filterCol(...
%     filename, ...
%     delim, ...
%     py.list({int32(1), int32(2)}), ...      % columns: chrom and pos
%     py.list({'eq','eq'}), ...              % operations: equality for both
%     py.list({1, 284844}), ...             % threshold values: 1 and 284844
%     skip, ...
%     bench, ...
%     py.None, ...                        % retrunCols = None: return full row
%     blocksize);
% result2 = char(result_py);
% disp('Result from filterCol (filtering chrom and pos):');
% disp(result2);

%% Test 3: Shutdown the Dask cluster
shutdownMsg = mod.bFileReaderDep.shutdownCluster();
disp('Cluster Shutdown Message:');
disp(char(shutdownMsg));

end