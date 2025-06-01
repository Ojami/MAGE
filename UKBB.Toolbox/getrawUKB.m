function UKB_df = getrawUKB(df, basket, toTable)
% This function is basically equivalent to UKBtraitParser but in simpler
% and faster form, and mainly created for getting single data field. It
% does not perform data pruning (e.g. removing missing variables, convert
% to double, etc.). All variables are originally stored as "string", so
% post-processing is not a goal of this function.
% 
% INPUTS:
%       df: UK Biobank(UKB) data-field number in UKB: e.g. 4289. If
%           interested in only specific array/instance specify it instead;
%           e.g.4260-2.0
%       toTable (OPTIONAL): a logical: to convert to table (true) or struct
%                           (false, default).
% OUTPUT:
%       UKB_df: UK Biobank data-field data, with missing values being
%       replaced by <missing>.

% IMPORTANT: MATLAB version must be >= R2016b
% Oveis Jamialahmadi, August 2019, GU.

warning("This function is deprecated as UKB-RAP replaced old basket-type files from UKBB!")

if nargin < 3
    toTable = false;
end
if nargin < 1
    fprintf('data-field cannot be empty!\n');
    help getrawUKB.m
    return
end

% Convert df(data-field) to readable format
df = strcat('x', df, '_');
df = regexprep(regexprep(df, '[.-]', '_'), '_$', '');

wd = fileparts(which('getrawUKB.m'));
targetDir = fullfile(wd, basket); % Target directory

% get EID
try
    UKB_eid = load(fullfile(targetDir, 'UKB_eid.mat'));
catch
    targetDir = fullfile(wd, "UKBFileParser_"+basket); % Target directory
    UKB_eid = load(fullfile(targetDir, 'UKB_eid.mat'));
end
UKB_df.eid = double(UKB_eid.UKB_eid);

variableMapper = load(fullfile(targetDir, 'variableMapper')); % Mapper file
variableMapper = variableMapper.variableMapper;

fileIdx  = unique(variableMapper(contains(variableMapper(:,2), df), 1));
UDI_codes = variableMapper(contains(variableMapper(:,2), df), 2);

if isempty(fileIdx)
    error('oops!query data-field does not exist in this data-basket!')
end

% create a loading bar 
progressBar = floor(linspace(1, numel(fileIdx)+1, 11));
fprintf('loading UKB data [           ]')

% load over chunk MAT files
for ii = 1:numel(fileIdx)
    UDI_idx = variableMapper(ismember(variableMapper(:,1), ...
        fileIdx{ii}), 2);
    UDI_idx = intersect(UDI_idx, UDI_codes);
    rawData = load(fullfile(targetDir, ...
        ['UKB_', fileIdx{ii}, '.mat']), UDI_idx{:});
    f_names = fieldnames(rawData);
    
    % pull variables together.
    for jj = 1:numel(f_names)
        rawData.(f_names{jj})(rawData.(f_names{jj}) == "") = missing;
        UKB_df.(f_names{jj}) = rawData.(f_names{jj});
    end
    clear rawData UDI_idx
    
    progressTxt = [repmat('=', 1, sum(ii >= progressBar)),'>',...
            repmat(' ', 1, 10-sum(ii >= progressBar))];
    fprintf(repmat('\b', 1, 12))
    fprintf('%s]', progressTxt)
end
fprintf('\n')

if toTable % this is slower but easier for data handling.
    UKB_df = struct2table(UKB_df);
end

end % END

