function bim2mat(files)
% converts bim files to mat files for easier use.
% Oveis Jamialahmadi, GU, Oct. 2020.

if nargin < 1 
    fprintf('ERROR: missing input bim files!\n')
    return
elseif any(cellfun(@isempty, regexp(files, '.bim$')))
    fprintf('ERROR:input not a bim file!\n')
    return
end

for i = 1:numel(files)
    
    fprintf('%d(of %d)-reading %s\n', i, numel(files), files{i})
    tic
    bim = datastore(files{i}, 'FileExtensions', '.bim',...
        'TextType', 'string', 'Type', 'tabulartext');
    bim.SelectedVariableNames = {'Var1', 'Var2', 'Var4', 'Var5', 'Var6'};
    bim = readall(bim);
    bim.Properties.VariableNames = {'chr', 'snp', 'pos', 'a1', 'a2'};
    bim = table2struct(bim, 'ToScalar', true);
    save(regexprep(files{i}, 'bim$' ,'mat'), '-struct', 'bim');
    toc
    
end

fprintf('Done!\n')
end % END