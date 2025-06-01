function a = struct2array(s)

    % see here: https://se.mathworks.com/matlabcentral/answers/1717910-was-struct2array-removed-from-matlab
    c = struct2cell(s);
    a = vertcat(c{:});

end