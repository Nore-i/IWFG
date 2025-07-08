% STK_IWFG returns bottom K hypervolume contributor from a list of points
%
% CALL: indices = stk_IWFG(Y, Y_REF, K)

function result = stk_IWFG (y, y_ref, k)

disp(y);
disp(y_ref);
disp(k);

% Missing or empty y_ref: will use [0 0 ... 0] as a reference point
if nargin < 2
    y_ref = [];
    k = 1;
end

if ~ isempty (y_ref)
    y_ref = double (y_ref);
    d = size (y_ref, 2);
    if (d == 0) || (~ isequal (size (y_ref), [1 d]))  % isrow CG#08
        stk_error ('y_ref should be a numeric row vector.', 'IncorrectSize');
    end
end

if nargin < 3
    k = 1;
else
    k = int32 (k);
end

% Pre-processing
if iscell (y),
    y = cellfun (@(z) wfg_preprocessing (z, y_ref), y, 'UniformOutput', false);
else % y is a matrix
    try
        y = double (y);
        assert (ndims (y) == 2);  %#ok<ISMAT> see CODING_GUDELINES
    catch
        stk_error (['y should either be a cell array or be (convertible ' ...
            'to) a numeric matrix'], 'InvalidArgument');
    end
    y = wfg_preprocessing (y, y_ref);
end

% Compute the hypervolume only
result = stk_IWFG_mex (y, k);

end % function


function y = wfg_preprocessing (y, y_ref)

y = double (y);

% Keep only non-dominated points, and remove duplicates
idx  = stk_paretofind(y);
mask = false(size(y,1),1);
mask(idx) = true;
y    = y(mask,:);
y    = unique(y,'rows','stable');

if isempty (y_ref)  % Use [0 0 ... 0] as a reference point
    
    % WFG convention: maximization problem
    y = - y;
    
else  % Reference point provided
    
    p = size (y, 2);
    p_ref = size (y_ref, 2);
    
    % Check the size of y
    if (p > p_ref)
        stk_error (['The number of columns the data matrix should not be ' ...
            'larger than the number of columns of y_ref'], 'InvalidArgument');
    end
    
    % WFG convention: maximization problem
    y = bsxfun (@minus, y_ref(1:p), y);
    
end

% Remove points that do not dominate the reference
b = any (y < 0, 2);
y(b, :) = [];

end % function



%!error hv = stk_dominatedhv ();
%!error hv = stk_dominatedhv (-y, 'incorrect ref type');
%!error hv = stk_dominatedhv (-y, [0 0]);
%!error hv = stk_dominatedhv (-y, [0 0 0 0 0], 0, 'too many input args');

%-------------------------------------------------------------------------------

%!shared y, hv0 % Example from README.TXT in WFG 1.10
%!
%! y = [ ...
%!     0.598 0.737 0.131 0.916 6.745; ...
%!     0.263 0.740 0.449 0.753 6.964; ...
%!     0.109 8.483 0.199 0.302 8.872 ];
%!
%! hv0 = 1.1452351120;

%!test
%! hv = stk_dominatedhv (-y);
%! assert (stk_isequal_tolrel (hv, hv0, 1e-10));

%!test
%! yy = stk_dataframe (- y);  % Check that @stk_dataframe inputs are accepted
%! hv = stk_dominatedhv (yy);
%! assert (stk_isequal_tolrel (hv, hv0, 1e-10));

%!test
%! hv = stk_dominatedhv (-y, [], 0);
%! assert (stk_isequal_tolrel (hv, hv0, 1e-10));

%!test
%! hv = stk_dominatedhv (-y, [0 0 0 0 0]);
%! assert (stk_isequal_tolrel (hv, hv0, 1e-10));

%!test
%! hv = stk_dominatedhv (1 - y, [1 1 1 1 1]);
%! assert (stk_isequal_tolrel (hv, hv0, 1e-10));

%!test
%! r = stk_dominatedhv (-y, [], 1);
%! hv = sum (r.sign .* prod (r.xmax - r.xmin, 2));
%! assert (stk_isequal_tolrel (hv, hv0, 1e-10));

%-------------------------------------------------------------------------------

%!shared y1, y2, y0, S, S1
%! y0 = [1.00 1.00];  % Reference point
%! y1 = [1.50 1.50];  % Above the reference point
%! y2 = [0.50 0.50];  % Below the reference point

%!assert (isequal (0.00, stk_dominatedhv (y1, y0)));
%!assert (isequal (0.25, stk_dominatedhv (y2, y0)));
%!assert (isequal (0.25, stk_dominatedhv ([y1; y2], y0)));
%!assert (isequal (0.25, stk_dominatedhv ([y2; y1; y2], y0)));

% Check decompositions:

%!test S = stk_dominatedhv (y1, y0, 1);    % empty decomposition
%!assert (isequal (size (S.xmin, 1), 0));

%!test S = stk_dominatedhv (y2, y0, 1);    % trivial decomposition
%!assert (isequal (S.sign, 1));
%!assert (isequal (S.xmin, y2));
%!assert (isequal (S.xmax, y0));

%!test S1 = stk_dominatedhv ([y2; y0], y0, 1);  % shoud be the same as before
%!assert (isequal (S1, S));

%!test S1 = stk_dominatedhv ([y2; y1], y0, 1);  % shoud be the same as before
%!assert (isequal (S1, S));

%!test S1 = stk_dominatedhv ([y2; y2], y0, 1);  % shoud be the same as before
%!assert (isequal (S1, S));

%-------------------------------------------------------------------------------

%!test
%! for d = 1:10,
%!    y = - 0.5 * ones (1, d);
%!    hv = stk_dominatedhv (y);
%!    assert (isequal (stk_dominatedhv (y), 0.5 ^ d));
%! end

%!test
%! for d = 1:10,
%!    y = - 0.5 * ones (1, d);
%!    r = stk_dominatedhv (y, [], 1);
%!    hv = sum (r.sign .* prod (r.xmax - r.xmin, 2));
%!    assert (isequal (stk_dominatedhv (y), 0.5 ^ d));
%! end

%-------------------------------------------------------------------------------

%!shared y, y_ref, dv, hv0
%! y1 = [0.25 0.75];
%! y2 = [0.50 0.50];
%! y3 = [0.75 0.25];
%!
%! y_ref = [1 1];
%!
%! y = {[], y1, y2, y3; [y1; y2], [y1; y3], [y2; y3], [y1; y2; y3]};
%!
%! dv = 0.25 ^ 2;  hv0 = [0 3 4 3; 5 5 5 6] * dv;

%!test
%! hv1 = stk_dominatedhv (y, y_ref);
%! assert (isequal (hv0, hv1));

%!test
%! r = stk_dominatedhv (y, y_ref, 1);
%!
%! % Check the first decomposition, which should be empty
%! assert (isempty (r(1).sign));
%! assert (isempty (r(1).xmin));
%! assert (isempty (r(1).xmax));
%!
%! % Check the other decompositions
%! for i = 2:6,
%!    hv2 = sum (r(i).sign .* prod (r(i).xmax - r(i).xmin, 2));
%!    assert (isequal (hv0(i), hv2));
%! end

%-------------------------------------------------------------------------------

%!test
%! y = (0.3:0.1:0.8)';
%! hv0 = 0.7;
%! hv1 = stk_dominatedhv (y, 1);
%! r = stk_dominatedhv (y, 1, true);
%! hv2 = sum (r.sign .* prod (r.xmax - r.xmin, 2));

%!test % four non-dominated points (hypervolume)
%! zr = [1 1];
%! zi = [0.2 0.8; 0.4 0.6; 0.6 0.4; 0.8 0.2]
%! P = perms (1:4);
%! for i = 1:24
%!     HV = stk_dominatedhv (zi(P(i, :), :), zr, 0);
%!     assert (stk_isequal_tolrel (HV, 0.4, 1e-15));
%! end

%!test % four non-dominated points (decomposition)
%! zr = [1 1];
%! zi = [0.2 0.8; 0.4 0.6; 0.6 0.4; 0.8 0.2]
%! P = perms (1:4);
%! for i = 1:24
%!     S = stk_dominatedhv (zi(P(i, :), :), zr, 1);
%!     HV = sum (S.sign .* prod (S.xmax - S.xmin, 2));
%!     assert (stk_isequal_tolrel (HV, 0.4, 1e-15));
%! end

%!test  % a case with 8 points and 5 objectives
%!      % http://sourceforge.net/p/kriging/tickets/33
%!
%! yr = [1.03 0.91 0.96 1.99 16.2];
%!
%! y = [ ...
%!     0.8180    0.5600    0.1264    1.0755    1.2462; ...
%!     0.8861    0.6928    0.2994    0.7228    0.9848; ...
%!     0.9021    0.8829    0.6060    0.1642    0.4282; ...
%!     0.9116    0.3097    0.8601    0.0468    0.2813; ...
%!     0.9306    0.1429    0.6688    0.1462    1.3661; ...
%!     0.9604    0.3406    0.4046    0.7239    1.8741; ...
%!     0.9648    0.7764    0.5199    0.4098    1.3436; ...
%!     0.9891    0.4518    0.7956    0.1164    1.2025];
%!
%! hv1 = stk_dominatedhv (y, yr, 0);
%!
%! S = stk_dominatedhv (y, yr, 1);
%! hv2 = sum (S.sign .* prod (S.xmax - S.xmin, 2));
%!
%! assert (isequal (size (S.sign), [87 1]));
%! assert (isequal (size (S.xmin), [87 5]));
%! assert (isequal (size (S.xmax), [87 5]));
%! assert (stk_isequal_tolrel (hv1, 1.538677420906463, 2 * eps));
%! assert (stk_isequal_tolrel (hv1, hv2, eps));

%!test % with random data
%! NREP = 5;
%! for p = 1:5
%!     for n = 1:10
%!         for i = 1:NREP
%!             % Draw random data
%!             y = rand (n, p);
%!             y = - y ./ (norm (y));
%!             % Compute hypervolume directly
%!             hv1 = stk_dominatedhv (y, [], 0);
%!             % Compute decomposition, then hypervolume
%!             R = stk_dominatedhv (y, [], 1);
%!             hv2 = sum (R.sign .* prod (R.xmax - R.xmin, 2));
%!             % Compare results
%!             assert (stk_isequal_tolabs (hv1, hv2, eps));
%!         end
%!     end
%! end
