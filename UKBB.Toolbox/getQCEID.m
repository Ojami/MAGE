function QCEID = getQCEID(qcparam, verbose, homedir)
% This function is part of UKBpheWAS to parse QC_EID.mat file under
% different scenarios
% qcparam: numeric for different scenarios:
%          1- s1_1, s2_1 and rest are merged.
%          2- s1_1, s2_2 and rest are merged.
%          3- s1_2, s2_1 and rest are merged.
%          4- s1_2, s2_2 and rest are merged.
% for details look at field .desc in QCfile
% 8/27/2019
% 8/4/2021: qcparam numeric was changed to string for better readability.
% 27/01/2022: a bug with verbose was fixed.
% 16/03/2022: new withdrawn samples were excluded.


if nargin < 1
    qcparam = 1;
    hlp = table(([-2, 0:6]).', ["All ethnicities"; ...
        "European & at most 10 putative third-degree relatives";...
        "European & maximum unrelated individuals"; ...
        "European & one member of each set of individuals with KING-estimated kinship coefficient >0.0442 is removed";...
        "white-British & maximum unrelatedness"; ...
        "white-British & one member of each set of individuals with KING-estimated kinship coefficient >0.0442 is removed";...
        "African & one member of each set of individuals with kinship coefficient >0.0442 is removed";...
        "Non-european & one member of each set of individuals with kinship coefficient >0.0442 is removed"],...
        'VariableNames', {'qcparam', 'meaning'});
    disp(hlp)
    disp('qcparam was set to 1 (AKA QC-MAX)!')
end

if nargin < 2
    verbose = true;
end

if nargin < 3
    homedir = fileparts(mfilename("fullpath"));
end

if ~isfile(fullfile(homedir,'QC_EID.mat'))
    error('getQCEID:cannot find QC_EID.mat in %s', homedir)
end

QCfile = load(fullfile(homedir,'QC_EID.mat'));
QCfile = QCfile.QC_EID;

switch qcparam
    case -2
        QCEID = QCfile.all;
        return
    case 0 % QC: used in mixed-models (BOLT and SAIGE)
        if verbose; fprintf('QC: European & at most 10 putative third-degree relatives\n'); end
        QCEID = setdiff(QCfile.european_naele, QCfile.excess_relatives);
    case 1 % QC-MAX:European & maximum unrelated individuals
        if verbose; fprintf('QC: European & maximum unrelated individuals\n'); end        
        QCEID = setdiff(QCfile.european_naele, union(QCfile.max_unrelated_neale, QCfile.excess_relatives));
    case 2 % European & one member of each set of individuals with KING-estimated kinship coefficient >0.0442 is removed
        if verbose; fprintf('QC: European & one member of each set of individuals with KING-estimated kinship coefficient >0.0442 is removed\n'); end
        QCEID = setdiff(QCfile.european_naele, QCfile.kinship_onepair);
    case 3 % white-British & at most 10 putative third-degree relatives
        if verbose; fprintf('QC: white-British & maximum unrelatedness\n'); end
        QCEID = setdiff(QCfile.white_british, QCfile.excess_relatives);
        QCEID = intersect(QCEID, QCfile.used_in_pca); % used.in.pca.calculation
    case 4 % white-British & one member of each set of individuals with KING-estimated kinship coefficient >0.0442 is removed
        if verbose; fprintf('QC: white-British & one member of each set of individuals with KING-estimated kinship coefficient >0.0442 is removed\n'); end
        QCEID = setdiff(QCfile.white_british, QCfile.kinship_onepair);
        QCEID = intersect(QCEID, QCfile.used_in_pca); % used.in.pca.calculation
    case 5 % Black 
        if verbose; fprintf('QC: African & one member of each set of individuals with kinship coefficient >0.0442 is removed\n'); end
        QCEID = setdiff(QCfile.black.eid, QCfile.black.kinship);
    case 6 % Non-europeans: setdiff of all participants and europeans
        if verbose; fprintf('QC: Non-european & one member of each set of individuals with kinship coefficient >0.0442 is removed\n'); end
        QCEID = setdiff(QCfile.noneuropean_neale.eid, QCfile.noneuropean_neale.kinship);
    otherwise
        error('Wrong input')
end

% apply exclusion criteria:
%   -withdrawn samples
%   -putative.sex.chromosome.aneuploidy
%   -sex mismatch
%   -Identified by the UK Biobank as outliers based on either genotyping
%    missingness rate or heterogeneity

withdrawn_fields = string(fieldnames(QCfile));
withdrawn_fields = withdrawn_fields(withdrawn_fields.lower.contains("withdrawn"));
exclude_fis = union(["sex_aneuploidy", "sex_mismatch", "ukb_outliers"], ...
   withdrawn_fields);

for k = 1:numel(exclude_fis)
    QCEID = setdiff(QCEID, QCfile.(exclude_fis(k))); 
end

end % END