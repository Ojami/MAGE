classdef pheno %< handle
% pheno calss implements methods for handling phenotype data (designed
% mainly for UK Biobank phenotypes, but generally can be applied to other
% cases) in a consistent manner.
% 
% Oveis Jamialahmadi, University of Gothenburg, December 2021.

   properties(Dependent)
       numeric
   end

   properties(GetAccess = public, SetAccess = private)
       eid
       value
   end

   properties
       term
       tag
   end
   
   % public ---------------------------------------------------------------
   methods
       function obj = pheno(eid, value)
           if nargin > 0
               if nargin < 2
                   error('pheno:MissingInput', 'both eid and value must be set')
               else
                   obj = checkSize(obj, eid, value);
               end
           end
       end
       
       function obj = fill(obj, eid, value)
           obj = checkSize(obj, eid, value);
       end
       
       function obj = keepEid(obj, eidN)
           obj = filterEid(obj, eidN, 'keep');
       end
       
       function obj = removeEid(obj, eidN)
           obj = filterEid(obj, eidN, 'remove');
       end
       
       function obj = filter(obj, cutoff, oper)
           arguments
               obj
               cutoff (1, 1) {mustBeNumeric}
               oper (1, 1) string {mustBeMember(oper, ["<=", "<", "=", ">", ">="])}
           end
           if obj.numeric
               operator.list = ["<=", "<", "=", ">", ">="];
               operator.fun = ["le", "lt", "eq", "gt", "ge"];
               oper = operator.fun(operator.list == oper);
               idx = feval(oper, obj.value, cutoff);
               obj = cleanPheno(obj, idx);
           else
               warning('method filter doesn''t work for non-numeric values!')
           end
       end
       
      function numeric = get.numeric(obj)
          % use class
           if all(isnumeric(obj.value))
               numeric = true;
           else
               numeric = false;
           end
      end
       
       function summary = getSummary(obj)
           if obj.numeric
               summary = valsummary(obj);
           else
               summary = catsummary(obj);
           end
       end
   end
   
   % private methods ------------------------------------------------------
   methods(Access = private)
       function s = valsummary(obj)
           s.mean = mean(obj.value);
           s.median = median(obj.value);
           s.iqr = iqr(obj.value);
       end
       
       function s = catsummary(obj)
           s.n = numel(obj.value);
       end

       function obj = checkSize(obj, eid, value)
           if numel(eid) ~= numel(value)
               error('pheno:SizeMismatch', 'numel eid must be equal to numel value!')
           elseif ~isvector(eid) || ~isvector(value)
               error('pheno:Vector', 'eid and values must be vectors!')
           else
               % check missingness: rmmissing
               if ~iscolumn(eid);   eid = eid.';     end
               if ~iscolumn(value); value = value.'; end
           end
           
           obj.eid = eid; obj.value = value;
       end
       
       function obj = filterEid(obj, eidN, C)
           idx = ismember(obj.eid, eidN);
           if strcmp(C, 'keep')
                idx = ~idx;
           end
           obj = cleanPheno(obj, idx);
       end
       
       function obj = cleanPheno(obj, idx)
           obj.eid(idx) = []; obj.value(idx) = [];
           if ~obj.numeric; obj.term(idx) = []; end
       end
       
   end
end