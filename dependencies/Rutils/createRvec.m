function vout = createRvec(v, n, str)

if nargin < 3, str = false; end
if nargin < 2, n = 10; end

% break a long vector into multiple lines to avoid wrapping 
step_size = min(numel(v), n);
idx = unique([1:step_size:(numel(v) + 1), numel(v) + 1]);

v = string(v);
if str
    ostr = 'c("';
    str = '","';
    endstr = '")';
    sep = '",';
    wrap = '"';
else
    ostr = "c(";
    str = ",";
    endstr = ")";
    sep = ",";
    wrap = "";
end

vout = strings(numel(idx) - 1, 1);
for k = 1:numel(idx) - 1
    tmp = join(v( idx(k):(idx(k+1) - 1) ), str);

    if numel(idx) <= 2
        vout(k) = ostr + tmp + endstr;
    else
        if k == 1
            vout(k) = ostr + tmp + sep;
        elseif k == (numel(idx) - 1)
            vout(k) = wrap + tmp + endstr;
        else
            vout(k) = wrap + tmp + sep;
        end
    end
end