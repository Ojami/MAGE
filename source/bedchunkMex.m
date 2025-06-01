function getBED = bedchunkMex(bedfile, nstart0, snplength, inSnpid, parflag)
if nargin < 5
    parflag = false;
end

if parflag
    getBED = cell(numel(inSnpid), 1);
    parfor ii = 1:numel(inSnpid)
        fid = fopen(bedfile, 'r');
        nstart = snplength*(inSnpid(ii)-1) + nstart0;
        nend = snplength*inSnpid(ii) + nstart0 - 1;
        fseek(fid, nstart-1, 'bof');  
        getbin = fread(fid,nend+1-nstart,'*int8');
        getbin = [bitget(getbin, 1), bitget(getbin, 2), bitget(getbin, 3), ...
            bitget(getbin, 4), bitget(getbin, 5), bitget(getbin, 6),...
            bitget(getbin, 7), bitget(getbin, 8)]';

        getBED{ii} = [getbin(1:2:end)', getbin(2:2:end)'];
        getBED{ii}(getbin(1:2:end) & ~getbin(2:2:end), :) = -2;
        fclose(fid);
    end
else

    fid = fopen(bedfile, 'r');
    getBED = cell(numel(inSnpid), 1);
    for ii = 1:numel(inSnpid)
        nstart = snplength*(inSnpid(ii)-1) + nstart0;
        nend = snplength*inSnpid(ii) + nstart0 - 1;
        fseek(fid, nstart-1, 'bof');  
        getbin = fread(fid,nend+1-nstart,'*int8');
        getbin = [bitget(getbin, 1), bitget(getbin, 2), bitget(getbin, 3), ...
            bitget(getbin, 4), bitget(getbin, 5), bitget(getbin, 6),...
            bitget(getbin, 7), bitget(getbin, 8)]';
        
        getBED{ii} = [getbin(1:2:end)', getbin(2:2:end)'];
        getBED{ii}(getbin(1:2:end) & ~getbin(2:2:end), :) = -2;
    end
    fclose(fid);
end

getBED = horzcat(getBED{:});

end % END