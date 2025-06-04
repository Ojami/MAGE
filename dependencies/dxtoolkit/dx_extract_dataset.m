function csv_files = dx_extract_dataset(df, opts)
% a wrapper for extract_dataset and tab-explorer (app) on UKB-RAP to
% extract phenotypes. Ideally, this function should be used in combination
% with phenoParser for follow-up disease endpoint definitions.
% For few phenotypes use extract_dataset, and for larger/complex (e.g.
% HESIN) use tab-explorer.
% 
% Oveis Jamialahmadi, University of Gothenburg, September 2024.

arguments
    df {mustBeA(df, "table")} % a dictionary from getUKBdictionary output for which fields should be fetched
    opts.instance_type {mustBeTextScalar} = "mem2_ssd1_v2_x16" % for table-exporter. see dx run --instance-type-help
    opts.list_instance_type (1,1) logical = false % if true only lists the available instances on UKB-RAP
    opts.dir {mustBeTextScalar} = fullfile(fileparts(which("phenoParser.m")), "UKBFileParser_RAP") % output dir, if does not exist, will create one.
    opts.updatefreq (1,1) double = 30 % update frequency for HESIN and DEATH tables
    opts.fetchOnly (1,1) logical = false % if true, only fetches 'df' without parsing (csv files will be written to opts.dir)

    % for UKB_HESIN_DEATH_2mat/UKBBasketParser
    opts.workers (1,1) double = 35 % number of workers.
end

% df = load("tdic.mat").tdic;

% Help for instant types
if opts.list_instance_type
    runbash("dx run --instance-type-help > dx_instance.txt");
    msg = readlines("dx_instance.txt");
    delete("dx_instance.txt")
    msg(~msg.startsWith("│"), :) = [];
    msg = split(msg, "│");
    msg(:, [1, end]) = [];
    msg = strip(msg);
    hdr = msg(1, :); msg(1, :) = [];
    msg = array2table(msg, VariableNames=hdr);
    msg = convertvars(msg, 2:width(msg), @double);
    disp(msg);
    return
end

dict = getUKBdictionary(ukbrap=true); % load UKB-RAP dictionary

if ~isfolder(opts.dir)
    mkdir(opts.dir);
end

% check variableMapper in opts.dir and update outdated (opts.updatefreq)
% fields
vmap_file = fullfile(opts.dir, "variableMapper.mat");

% check if it's empty --> delete it
if isfile(vmap_file)
    vmap = load(vmap_file);
    if isempty(vmap.datev)
        delete(vmap_file);
    end
end

if isfile(vmap_file) % otherwise it's the first try to fetch data
    vmap = load(vmap_file);
    datev = vmap.datev;
    variableMapper = vmap.variableMapper;
    vtab = array2table(variableMapper, VariableNames=["chunk", "name"]);
    vtab.date = datev;

    % outdated fields
    vtab.update = days(datetime("today") - vtab.date) > opts.updatefreq;

    % remove intersection of vtab updated fields and df fields
    df(ismember(df.name, vtab.name(~vtab.update)), :) = [];

    if isempty(df) % everything is upToDate
        return
    end

    % keep outdated variables and update if they're in df
    vtab(~vtab.update & ~ismember(vtab.name, df.name), :) = [];
    if ~isempty(vtab)
        % delete outdated fields from variableMapper. Note that, it's more
        % efficient to work with matfile, so here by 'delete' we mean
        % replace with an empty vector, though the variable (emptied) is
        % still present in that mat file. But since we update
        % variableMapper, phenoParser.m function does not use that variable
        % anymore.
        
        uchunk = unique(vtab.chunk);
        for k = 1:numel(uchunk)
            idx = vtab.chunk == uchunk(k);
            vtmp = vtab(idx, :);

            chunk_file = matfile(fullfile(opts.dir, "UKB_" + uchunk(k) + ".mat"), ...
                Writable=true);
            for j = 1:height(vtmp)
                chunk_file.(vtmp.name(j)) = [];
            end
            clear("chunk_file")

        end

        % update variableMapper file
        idx = ismember(variableMapper(:, 1), vtab.chunk) & ismember(variableMapper(:, 2), vtab.name);
        variableMapper(idx, :) = [];
        datev(idx, :) = [];
        save(vmap_file, "variableMapper", "datev")

    end
    

end

% if any of entities contains "hesin" or "death", jobs will be sent to
% table-exporter if  will be saved to UKB_HESIN_DEATH_2mat folder.
% Then UKB_HESIN_DEATH_2mat.m function is
% called internally to update the HESIN and DEATH mat files if data in the
% HESIN_DEATH folder is outdate (above opts.updatefreq days)
opts.gethesin = false;
hesin_death_idx = ismember(df.entity, ["hesin", "death", "death_cause", "hesin_diag"]);
if any(hesin_death_idx)
    df(hesin_death_idx, :) = []; % don't fetch these data. Will handle differently

    opts.hesindir = fileparts(which("UKB_HESIN_DEATH_2mat.m"));
    vfile = fullfile(opts.hesindir, "variableMapper.mat");
    if ~isfile(vfile)
        opts.gethesin = true;
    else
        vdate_old = load(vfile).datev; % last update of HESIN/DEATH data
        vdate_now = datetime("today");
        if days(vdate_now - vdate_old) > opts.updatefreq
            opts.gethesin = true;
        end
    end

    if opts.gethesin 
        opts.hesin = ["death" , "death_cause", "hesin", "hesin_diag"];
    end

end

if isempty(df) && ~ opts.gethesin % nothing to update or fetch
    return
end

if ~opts.dir.endsWith(filesep), opts.dir = opts.dir + filesep; end

% find project ID
[~, project_id] = runbash("dx env | grep project- | awk -F '\t' '{print $2}'");
project_id = project_id(end);
[~, dataset_name] = runbash("dx ls " + project_id, wait=true);
dataset_name(~dataset_name.endsWith(".dataset")) = [];
dataset_name = natsort(dataset_name);
dataset_name = dataset_name(end); % latest dataset

[~, dataset_id] = runbash("dx describe " + project_id + ":" + dataset_name);
dataset_id(~dataset_id.startsWith("ID")) = [];
dataset_id = dataset_id.extractAfter("ID").strip;


df = df(:, colnames(dict));
dict_eid = dict(dict.name == "eid", :);
uentity = unique(df.entity);
jobs = struct;
for k = 1:numel(uentity)
    
    % if entity is anything but participant, add the full list of fields
    if uentity(k) ~= "participant"
        tmp = dict(dict.entity == uentity(k), :);
    else
        tmp = df(df.entity == uentity(k), :);
    end
    
    % add eid if necessary
    if ~any(tmp.name == "eid")
        tmp2 = dict_eid(dict_eid.entity == uentity(k), :);
        ds = [tmp; tmp2];
    else
        ds = tmp;
    end

    ds.id = ds.entity + "." + ds.name;
    
    % try extract_dataset (locally) if entity is participant or olink and
    % has either instances or arrays
    if uentity(k).startsWith(["olink_", "participant"])
        
        % exclude those fields without arrays and instances (e.g.
        % hospitalized data with complex structure)
        idx_nan = ismissing(ds.instance) & ismissing(ds.array) & ismissing(ds.primary_key_type); % 'primary_key_type' for eid
        
        if any(~idx_nan) % is there anything to be run using extract_dataset?

            ds_extract = ds(~idx_nan, :);
            if all(ds_extract.FieldID == "eid") % don't fetch eid alone
                ds_export = ds;
            else
                
                ds_extract = ds_extract_chunk(ds_extract);

                tt = tic;
                fprintf("fetching %d chunks locally...\n", numel(ds_extract))
                
                failed_ds_extract = cell(numel(ds_extract), 1);
                for j = 1:numel(ds_extract)
                    tmp_name = getRandomName("tmp.csv", 10);
                    if isfile(fullfile(opts.dir, tmp_name))
                        tmp_name = getRandomName("tmp.csv", 12); % generate a new name
                    end
        
                    cmd = join(ds_extract{j}.id, ",");
                    cmd = "dx extract_dataset " + project_id + ":" + dataset_id + ...
                        " --delimiter ',' --output '" + makeWSLpath(opts.dir) + ...
                        tmp_name + "' " + " --fields " + cmd;
                    [~, msg] = runbash(cmd, verbose=false, parallel=false, wait=true);
        
                    if numel(msg) > 1 && any(msg ~= "") % failed: try table-exporter 
                        fprintf("\tchunk %d failed (will try with table-exporter)\n", j)
                        failed_ds_extract{j} = ds_extract{j};
                    else
                        fprintf("\tchunk %d was fetched successfully.\n", j)
                    end

                end % loop over chunks of ds_extract (local feth)

                tt = toc(tt);
                tt = duration(0, 0, tt);
                fprintf('done in %s (hr:min:sec)\n', tt)

                % check failed extracts
                failed_ds_extract(cellfun(@isempty, failed_ds_extract)) = [];
                failed_ds_extract = vertcat(failed_ds_extract{:});
                if ~isempty(failed_ds_extract) % keep only one eid record
                    eid_idx = find(failed_ds_extract.name == "eid");
                    failed_eid = failed_ds_extract(eid_idx(1), :);
                    failed_ds_extract(eid_idx, :) = [];
                    failed_ds_extract_copy = [failed_eid; failed_ds_extract];

                    ds_export = [ds(idx_nan, :); failed_ds_extract_copy];
                else
                    ds_export = ds(idx_nan, :);
                end
                
            end

        end
        
        if ~isempty(ds_export) % is there anything to be sent to table-exporter? 
            jobs.tab{k} = ds_export;
        end

    else % send to table-exporter
        jobs.tab{k} = ds;
    end

    clear ds tmp
end

% set output folder for table-exporter
[~, files_main_path] = runbash("dx ls "+ project_id, wait=true);
if ~any(files_main_path == "table_exporter/")
    runbash("dx mkdir " + project_id + ":table_exporter");
end
opts.destination = project_id + ":table_exporter";


% send jobs to table-exporter
ct = 1;
if isfield(jobs, "tab")
    miss_idx = cellfun(@isempty, jobs.tab);
    jobs.tab(miss_idx) = [];

    if ~isempty(jobs.tab) % anything to fetch?
        for k = 1:numel(jobs.tab)
            df_in = jobs.tab{k};

            if df_in.entity(1) == "participant"
                df_in(df_in.name == "eid", :) = [];
    
                field_names = " -ifield_names=eid ";
                for j = 1:height(df_in)
                    field_names = field_names + "-ifield_names=" + df_in.name(j) + " ";
                end
    
                cmd = "dx run table-exporter" + ...
                    " -idataset_or_cohort_or_dashboard=" + dataset_id + ...
                    field_names + ...
                    " -ientity=" + df_in.entity(1) + ...
                    " -ioutput=tab" + k + ...
                    " -ioutput_format=CSV -iheader_style=FIELD-NAME" + ...
                    " -icoding_option=RAW" + ...
                    " --destination " + opts.destination + ...
                    " --instance-type " + opts.instance_type + ...
                    " --brief -y --ignore-reuse --tag tab" + k;
                [~, opts.jobid(ct)] = runbash(cmd, verbose=true); % send jobs
                ct = ct + 1;
    
                % df.tag = df.entity + "_" + df.name;
                % for j = 1:height(df) % different jobs per each complex df
                %     cmd = "dx run table-exporter" + ...
                %         " -idataset_or_cohort_or_dashboard=" + dataset_id + ...
                %         field_names + ...
                %         " -ientity=" + df.entity(j) + ...
                %         " -ioutput=" + df.tag(j) + ...
                %         " -ioutput_format=CSV -iheader_style=FIELD-NAME" + ...
                %         " -icoding_option=RAW" + ...
                %         " --destination " + opts.destination + ...
                %         " --instance-type " + opts.instance_type + ...
                %         " --brief -y --ignore-reuse --tag " + df.tag(j);
                %     [~, opts.jobid(ct)] = runbash(cmd, verbose=true); % send jobs
                %     ct = ct + 1;
                % end

            else % other entities such as omop, gp, ...: fetch whole entity table
                cmd = "dx run table-exporter" + ...
                    " -idataset_or_cohort_or_dashboard=" + dataset_id + ...
                    " -ientity=" + df_in.entity(1) + ...
                    " -ioutput=" + df_in.entity(1) + ...
                    " -ioutput_format=CSV -iheader_style=UKB-FORMAT" + ...
                    " -icoding_option=RAW" + ...
                    " --destination " + opts.destination + ...
                    " --instance-type " + opts.instance_type + ...
                    " --brief -y --ignore-reuse --tag " + df_in.entity(1);
                [~, opts.jobid(ct)] = runbash(cmd, verbose=true); % send jobs
                ct = ct + 1;
            end
        end
    end

end

% check if hesin tables should be fetched
if opts.gethesin
    for k = 1:numel(opts.hesin)
        cmd = "dx run table-exporter" + ...
            " -idataset_or_cohort_or_dashboard=" + dataset_id + ...
            " -ientity=" + opts.hesin(k) + ...
            " -ioutput=" + opts.hesin(k) + ...
            " -ioutput_format=CSV -iheader_style=UKB-FORMAT" + ...
            " -icoding_option=RAW" + ...
            " --destination " + opts.destination + ...
            " --instance-type " + opts.instance_type + ...
            " --brief -y --ignore-reuse --tag " + opts.hesin(k);
        [~, opts.jobid(ct)] = runbash(cmd, verbose=true); % send jobs
        ct = ct + 1;
    end
end

if isfield(opts, "jobid") && ~isempty(opts.jobid)

    tt = tic;
    fprintf("waiting for %d jobs on DNANexus...\n", numel(opts.jobid))
    allRunning = true;
    [total_failed_jobs, total_done_jobs] = deal(0); % keep track of failed/done jobs
    while allRunning

        opts.jobstat = cell(numel(opts.jobid), 1);
        for k = 1:numel(opts.jobid)
            [~, msg] = runbash("dx describe " + opts.jobid(k) + " --color off");
            msg(~msg.startsWith(["ID", "State", "Tags", "Output"])) = [];
            msg(msg.startsWith("Output folder")) = [];

            % extract file id
            fidx = msg.startsWith("Output");
            file_id = msg(fidx).extractBetween("[", "]").strip;
            if isempty(file_id)
                file_id = "-";
            end
            msg(fidx, :) = [];
            msg = msg.split';
            opts.jobstat{k} = array2table(msg(2, :), VariableNames=msg(1, :));
            opts.jobstat{k}.Output = file_id;
        end

        chekstats = vertcat(opts.jobstat{:});

        % save job details for download
        if ~isfield(opts, "job_tags")
            opts.job_tags = chekstats;
            opts.job_tags.done(:) = false;
        end

        done_jobs = chekstats.State.contains("done");
        if any(done_jobs)
            total_done_jobs = total_done_jobs + sum(done_jobs);
            fprintf("\t%d jobs are done so far.\n", total_done_jobs)
        end

        failed_jobs = chekstats.State.contains("failed");
        if any(failed_jobs)
            total_failed_jobs = total_failed_jobs + sum(failed_jobs);
            fprintf("\t%d jobs are failed so far.\n", total_failed_jobs)
        end

        skip_ids = chekstats.ID(failed_jobs | done_jobs);
        if ~isempty(skip_ids)
            opts.job_tags.done(ismember(opts.job_tags.ID, chekstats.ID(done_jobs))) = true;
            opts.jobid(ismember(opts.jobid, skip_ids)) = []; % don't check again

            if ~isempty(opts.jobid)
                fprintf("\tstill waiting for %d jobs.\n", numel(opts.jobid))
            end
        end

        if isempty(opts.jobid) % all jobs are finished (failed or done)
            allRunning = false;
        end

    end

    tt = toc(tt);
    tt = duration(0, 0, tt);
    fprintf('done in %s (hr:min:sec)\n', tt)

    % download and remove files from platform
    done_jobs = opts.job_tags(opts.job_tags.done, :);
    if ~isempty(done_jobs)
        
        % get file-ids
        jobs_tab = cell(height(done_jobs), 1);
        for k = 1:height(done_jobs)
            [~, msg] = runbash("dx describe " + done_jobs.ID(k) + " --color off");
            msg(~msg.startsWith(["ID", "State", "Tags", "Output"])) = [];
            msg(msg.startsWith("Output folder")) = [];

            % extract file id
            fidx = msg.startsWith("Output");
            file_id = msg(fidx).extractBetween("[", "]").strip;
            if isempty(file_id)
                file_id = "-";
            end
            msg(fidx, :) = [];
            msg = msg.split';
            jobs_tab{k} = array2table(msg(2, :), VariableNames=msg(1, :));
            jobs_tab{k}.Output = file_id;
        end

        jobs_tab = vertcat(jobs_tab{:});

        tt = tic;
        fprintf("downloading data from DNANexus...")

        for k = 1:height(done_jobs)
            cmd = "dx download " + jobs_tab.Output(k) + ...
                " --output " + makeWSLpath(opts.dir) + " --overwrite";
            runbash(cmd, verbose=true);

            cmd = "dx rm " + jobs_tab.Output(k);
            runbash(cmd, verbose=true);
        end

        tt = toc(tt);
        tt = duration(0, 0, tt);
        fprintf('\b\b done in %s (hr:min:sec)\n', tt)
    end

end

if opts.fetchOnly
    csv_files = getfilenames(opts.dir, "csv").csv;
    if ~isempty(csv_files)
        disp("following files have been fetched:")
        fprintf("\t%s\n", csv_files)
        csv_files = fullfile(opts.dir, csv_files);
    else
        csv_files = [];
    end
    return
else

    % call UKBBasketParser for remaining files
    % find all csv files in the opts.dir
    csv_files = getfilenames(opts.dir, "csv", "fullpath", false).csv;
end

% remove hesin files
csv_files(csv_files.endsWith(["death", "death_cause", "hesin", "hesin_diag"] + ".csv")) = [];

% remove other non "participant" datasets
idx = csv_files.endsWith(df.entity + ".csv");
csv_entity = csv_files(idx, :);
csv_files(idx, :) = [];

if ~isempty(csv_files)
    tt = tic;
    fprintf("Parsing %d participant datasets...", numel(csv_files))
    
    csv_files = fullfile(opts.dir, csv_files);
    for k = 1:numel(csv_files)
        UKBBasketParser(csv_files(k), ukbrap=true, threads=opts.workers, ...
            path=opts.dir, verbose=false)
        delete(csv_files(k))
    end

    tt = toc(tt);
    tt = duration(0, 0, tt);
    fprintf('\b\b done in %s (hr:min:sec)\n', tt)
end

% parse non-"participant" datasets
if ~isempty(csv_entity)
    tt = tic;
    fprintf("Parsing %d entity datasets...", numel(csv_entity))

    csv_entity = fullfile(opts.dir, csv_entity);
    for k = 1:numel(csv_entity)
        UKBBasketParser(csv_entity(k), ukbrap=true, threads=opts.workers, ...
            path=opts.dir, verbose=false, entity=true)
        delete(csv_entity(k))
    end

    tt = toc(tt);
    tt = duration(0, 0, tt);
    fprintf('\b\b done in %s (hr:min:sec)\n', tt)
end

% call UKB_HESIN_DEATH_2mat if necessary
if opts.gethesin
    tt = tic;
    fprintf("Parsing HESIN and DEATH datasets.\n")
    UKB_HESIN_DEATH_2mat(death=fullfile(opts.dir, "death.csv"), ...
        death_cause=fullfile(opts.dir, "death_cause.csv"), ...
        hesin=fullfile(opts.dir, "hesin.csv"), ...
        hesin_diag=fullfile(opts.dir, "hesin_diag.csv"), ...
        dir=opts.hesindir, workers=opts.workers)

    % delete residual files
    delete(fullfile(opts.dir, "death.csv"));
    delete(fullfile(opts.dir, "death_cause.csv"));
    delete(fullfile(opts.dir, "hesin.csv"));
    delete(fullfile(opts.dir, "hesin_diag.csv"));

    tt = toc(tt);
    tt = duration(0, 0, tt);
    fprintf('Done in %s (hr:min:sec)\n', tt)
end

end % END

%% subfunctions ===========================================================
function ds = ds_extract_chunk(df)
% split input dfs into chunks of dictionary tables to be fetched and avoid
% runtime errors.
cs = 15; % chunk size

ds = cell({});
if height(df) <= cs
    ds{1} = df;
    return
end

% get eid
idx = df.name == "eid";
eid_df = df(idx, :); df(idx, :) = [];

ch = unique([1:cs:height(df), height(df)+1]);
for k = 1:(numel(ch)-1)
    idx_range = ch(k):(ch(k+1)-1);
    ds{k} = [eid_df; df(idx_range, :)];
end

end % END
