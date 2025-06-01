function [ingwas, chr] = snptag(ingwas, bimfiles, verbose)
if nargin < 3
    verbose = false;
end

snpcol = find(contains(ingwas.Properties.VariableNames, {'SNP', 'snp', 'annotate', 'variant'}));
snpcol = snpcol(1); % there should be only one column for variant ids, right?

allelecols = find(ismember(ingwas.Properties.VariableNames, {'allele0', 'allele1'}));

bedext = endsWith(bimfiles, '.bed');
bimfiles(bedext) = regexprep(bimfiles(bedext), '.bed$', '.bim');
bimfiles(~bedext) = bimfiles(~bedext) + ".bim";

if ~isstring(bimfiles)
    bimfiles = string(bimfiles);
end

% get chromosomes in bim files
chr = strings(numel(bimfiles), 1);
for i = 1:numel(bimfiles)
    fid = fopen(bimfiles(i), 'r');
    bimLine = string(split(fgetl(fid)));
    chr(i) = bimLine(1);
    fclose(fid);
end

% rawBIM = [];
for ii = 1:numel(bimfiles)
    
    if verbose
        fprintf('%d(of %d)-checking/modifying variant IDs...\n', ii, numel(bimfiles))
    end

    % work only with variants of this bimfile 
    chridx = string(ingwas.chr) == chr(ii);
    gwas = ingwas(chridx, [snpcol, allelecols]);
    snp = find(contains(gwas.Properties.VariableNames, {'SNP', 'snp', 'annotate', 'variant'}));

    % read bim file
    getBIMraw = datastore(bimfiles(ii), 'Type', 'tabulartext',...
        'FileExtensions', '.bim', 'TextType', 'string');
    getBIMraw.SelectedVariableNames = {'Var2'};
    getBIMraw = readall(getBIMraw);
    getBIMraw = getBIMraw.Var2;
%     rawBIM = [rawBIM; getBIMraw];
    
    % identify tagged variants
    f_tagged = contains(getBIMraw, '_') & startsWith(getBIMraw, "rs");
    if ~any(f_tagged), continue; end %@08AUG2024
    getBIM = getBIMraw(f_tagged);
    
    % new approach --------------------------------------------------------
    % @12/01/2022: based on new modifications to bimtag function, following
    % steps are not needed anymore, and tags can be easily created (given
    % that the original gwas summary stat file was not modified).

    % different snps, either 1-gwas has fewer snps or 2-different names(tags)

    %@08AUG2024: a bug was fixed
    gwas.id1 = gwas.(snp) + "_" + gwas.allele0 + "_" + gwas.allele1;
    gwas.id2 = gwas.(snp) + "_" + gwas.allele1 + "_" + gwas.allele0;
    [idx1, idx2] = ismember(gwas.id1, getBIM);
    if any(idx1), gwas.(snp)(idx1) = getBIM(idx2(idx1)); end
    [idx1, idx2] = ismember(gwas.id2, getBIM);
    if any(idx1), gwas.(snp)(idx1) = getBIM(idx2(idx1)); end
    assert(all(ismember(gwas.(snp), getBIMraw)), "snptag: not all GWAS variants found in the reference panel!")
    ingwas.(snpcol)(chridx) = gwas.(snp);

    % dupsnps = ismember(gwas.(snp), duplicates(gwas.(snp)));
    % mixedtags = gwas.(snp)(dupsnps) + "_" + gwas.allele0(dupsnps) + "_" + gwas.allele1(dupsnps);
    % getBIM = setdiff(getBIM, gwas.(snp)); % tagged variants
    % if all(ismember(mixedtags, getBIM)) % all modified snp names could be found in mixedtags
    %     allfound = true;
    % else % trye allele1 + allele0
    %     mixedtags = gwas.(snp)(dupsnps) + "_" + gwas.allele1(dupsnps) + "_" + gwas.allele0(dupsnps);
    %     if all(ismember(mixedtags, getBIM)) % same as above
    %         allfound = true;
    %     else
    %         allfound = false; % go with the old way
    %     end
    % end
    % 
    % if allfound
    %     gwas.(snp)(dupsnps) = mixedtags;
    %     if ~any(ismember(gwas.(snp), getBIMraw))
    %         error('snptag: something went wrong!')
    %     end
    %     ingwas.(snpcol)(chridx) = gwas.(snp);
    %     continue
    % else
    %     error('snptag:should be fixed')
    % end
    
    % % old approach --------------------------------------------------------
    % % if above approach failed (old bim files) keep on with the old inefficient way
    % getBIM = setdiff(getBIM, gwas.(snp), 'stable'); % tagged variants
    % getBIM = regexp(getBIM, '_', 'split', 'once');
    % if any(cellfun(@isempty, getBIM))
    %     error('cannot find any tag in mismatched SNPs: %s', bimfiles(ii))
    % end
    % getBIM = vertcat(getBIM{:});
    % f_matched = ismember(gwas.(snp), getBIM(:, 1)); % matched variants
    % f_matched = find(f_matched);
    % checksnp = string(gwas.(snp)(f_matched));
    % checksnpU = unique(checksnp, 'stable');
    % 
    % if isempty(checksnp) % No need to modify modSNP
    %     continue
    % end
    % 
    % checkSNP = gwas.(snp);
    % for jj = 1:numel(checksnpU)
    %     idx2 = find(getBIM(:, 1) == checksnpU(jj));
    %     idx1 = f_matched(checksnp == checksnpU(jj));          
    %     tag =  gwas.allele1(idx1) + "_" + gwas.allele0(idx1);
    % 
    %     % there are exceptions like rs548234096 that cannot be identified
    %     % using the following condition, due to different MAF
    %     % distributions in qcmax (50 k subsample) and qc. So, before adding
    %     % numeric tags first check flipped versions.
    %     tag_flip = union(tag, gwas.allele0(idx1) + "_" + gwas.allele1(idx1));
    %     if numel(tag_flip) == 2 && numel(idx1) > 1 && isscalar(unique(tag))
    %         [f1, f2] = ismember(tag_flip, getBIM(idx2, 2));
    %         f2(f2<1) = [];
    %     else
    %         f1 = 0;
    %     end
    %     if any(~f1) % above condition failed
    %         if numel(tag) > numel(unique(tag))
    %             tag = tag + ["";(1:numel(tag)-1)'];
    %         end
    %         [f1, f2] = ismember(tag, getBIM(idx2, 2));
    %         if any(~f1) % missmatch found
    %             % possible reason: marginal MAF in reference BIM. Flip
    %             % and check again
    %             tag_flip =  gwas.allele0(idx1) + "_" + gwas.allele1(idx1); 
    %             if numel(tag_flip) > numel(unique(tag_flip))
    %                 tag_flip = tag_flip + (1:numel(tag_flip))';
    %             end
    %             [f1_f, f2_f] = ismember(tag_flip, getBIM(idx2, 2));
    % 
    %             if isempty(f2_f) || all(~f2_f)
    %                 error('cannot match tag:%s!\n', bimfiles(ii))
    %             end
    % 
    %             f1 = f1 | f1_f;
    %             f2 = f2 + f2_f;
    %         end 
    %         f2(f2<1) = [];
    %     end
    %     checkSNP(idx1(f1)) = cellstr(getBIM(idx2(f2), 1) + "_" + getBIM(idx2(f2), 2));
    % end
    % 
    % modSNP.(snp) = checkSNP; % incorrect! should be changed 
end

% if any(~ismember(modSNP.(snp), rawBIM))
%     error('Something went wrong!')
% end
end


%% ========================================================================
% function modSNP = snptag(modSNP, bimfiles, chrom)
% 
% if ~isstring(bimfiles)
%     bimfiles = string(bimfiles);
% end
% 
% bed_chrom = regexp(bimfiles, '\d+$', 'match');
% bed_chrom = double(string(bed_chrom));
% if ~isempty(bed_chrom) % otherwise: imp.lof
%     bimfiles = bimfiles(ismember(bed_chrom, chrom));
% end
% 
% for ii = 1:numel(bimfiles)
%     % read bim file
%     getBIMraw = datastore(bimfiles(ii)+".bim", 'Type', 'tabulartext',...
%         'FileExtensions', '.bim', 'TextType', 'string');
%     getBIMraw.SelectedVariableNames = {'Var2'};
%     getBIMraw = readall(getBIMraw);
%     getBIMraw = getBIMraw.Var2;
%     
%     % identify tagged variants
%     f_tagged = contains(getBIMraw, '_');
%     getBIM = getBIMraw(f_tagged);
%     getBIM = setdiff(getBIM, modSNP.SNP, 'stable'); % tagged variants
%     getBIM = regexp(getBIM, '_', 'split', 'once');
%     if any(cellfun(@isempty, getBIM))
%         error('cannot find any tag in mismatched SNPs: %s', bimfiles(ii))
%     end
%     getBIM = vertcat(getBIM{:});
%     f_matched = ismember(modSNP.SNP, getBIM(:, 1)); % matched variants
%     f_matched = find(f_matched);
%     checksnp = string(modSNP.SNP(f_matched));
%     checksnpU = unique(checksnp, 'stable');
%     
%     if isempty(checksnp) % No need to modify modSNP
%         continue
%     end
%     
%     checkSNP = modSNP.SNP;
%     for jj = 1:numel(checksnpU)
%         idx2 = find(getBIM(:, 1) == checksnpU(jj));
%         idx1 = f_matched(checksnp == checksnpU(jj));          
%         tag =  modSNP.A2(idx1) + "_" + modSNP.A1(idx1);
%         
%         % there are exceptions like rs548234096 that cannot be identified
%         % using the following condition, due to different MAF
%         % distributions in qcmax (50 k subsample) and qc. So, before adding
%         % numeric tags first check flipped versions.
%         tag_flip = union(tag, modSNP.A1(idx1) + "_" + modSNP.A2(idx1));
%         if numel(tag_flip) == 2 && numel(idx1) > 1 && numel(unique(tag)) == 1
%             [f1, f2] = ismember(tag_flip, getBIM(idx2, 2));
%             f2(f2<1) = [];
%         else
%             f1 = 0;
%         end
%         if any(~f1) % above condition failed
%             if numel(tag) > numel(unique(tag))
%                 tag = tag + ["";(1:numel(tag)-1)'];
%             end
%             [f1, f2] = ismember(tag, getBIM(idx2, 2));
%             if any(~f1) % missmatch found
%                 % possible reason: marginal MAF in reference BIM. Flip
%                 % and check again
%                 tag_flip =  modSNP.A1(idx1) + "_" + modSNP.A2(idx1); 
%                 if numel(tag_flip) > numel(unique(tag_flip))
%                     tag_flip = tag_flip + (1:numel(tag_flip))';
%                 end
%                 [f1_f, f2_f] = ismember(tag_flip, getBIM(idx2, 2));
% 
%                 if isempty(f2_f) || all(~f2_f)
%                     error('cannot match tag:%s!\n', bimfiles(ii))
%                 end
% 
%                 f1 = f1 | f1_f;
%                 f2 = f2 + f2_f;
%             end 
%             f2(f2<1) = [];
%         end
%         checkSNP(idx1(f1)) = cellstr(getBIM(idx2(f2), 1) + "_" + getBIM(idx2(f2), 2));
%     end
%     modSNP.SNP = checkSNP;
% end
% 
% if any(~ismember(modSNP.SNP, getBIMraw))
%     error('Something went wrong!')
% end
% end

%%
% for jj = 1:numel(checksnpU)
%     idx1 = find(getBIM(:, 1) == checksnpU(jj));
%     idx2 = f_matched(checksnp == checksnpU(jj));          
%     al1 = modSNP.A1(idx2); al2 = modSNP.A2(idx2);
% 
%     [f_curr1, f_ref1] = ismember(al1, getBIM(idx1, 2));
%     [f_curr2, f_ref2] = ismember(al2, getBIM(idx1, 2));
%     f_ref1(f_ref1<1) = []; f_ref2(f_ref2<1) = [];
% 
%     if ~all(f_curr1) && ~all(f_curr2) % cannot match
%         al1 = string(al1); al2 = string(al2);
%         al1_len = arrayfun(@(x)length(x{:}),al1);
%         al2_len = arrayfun(@(x)length(x{:}),al2);
%         if numel(unique([al1_len; al2_len])) < 2 %e.g. C/G and G/T
%             al_shared = intersect(al1,al2);
%             al1(ismember(al1, al_shared)) = "";
%             al2(ismember(al2, al_shared)) = "";
%             al1 = al1 + al2; 
%             al2 = repmat(al_shared, numel(al1), 1);
%         else % e.g. GC/GCC and G/GC
%             if isempty(setdiff(al1, al2)) % e.g. A/AT and AT/A
%                 al1 = string(1:numel(al1))';
%             else
%                 al12_diff = al1_len - al2_len;
%                 al12_diff(al12_diff >= 0) = 1; al12_diff(al12_diff < 0) = 0;
%                 al12_diff = logical(al12_diff);
%                 flip_ch = al1(al12_diff);
%                 al1(al12_diff) = al2(al12_diff); al2(al12_diff) = flip_ch;
%                 al1 = al1 + "_" + al2;
%             end
%         end
%         [f_curr1, f_ref1] = ismember(al1, getBIM(idx1, 2));
%         [f_curr2, f_ref2] = ismember(al2, getBIM(idx1, 2));
%         f_ref1(f_ref1<1) = []; f_ref2(f_ref2<1) = [];
%     end
%     idx1 = idx1([f_ref1; f_ref2]);
%     idx2 = idx2(f_curr1 | f_curr2);
% 
%     if ~all(strcmp(al1,getBIM(idx1, 2))) && ~all(strcmp(al2,getBIM(idx1, 2)))
%         fprintf('ERROR: cannot match to ref bim file!\n')
%         modSNP = [];
%         return
%     end
%     modSNP.SNP(idx2) = cellstr(getBIM(idx1, 1) + "_" + getBIM(idx1, 2));
% end