function [status, results] = runAll(varargin)
% Runs all tests in the 'toupee' repo
% Gathers and runs all tests in the recursive contents of this directory.
% All input arguments should be name-value pairs.
% 
% Inputs:
% -------
% 'ignore' : cell array of chars
%   The folder directories of tests to *not* include.      
% 
% Outputs:
% --------
% 'status' : logical
%   0 if all the tests that were run passed, 1 otherwise.
% 'results' : TestResult array
%   Full details on all the tests that were run.
%
% Examples:
% ---------
% 1) Run all the tests in the repo:
%   [status, failures] = toupee.tests.runAll();
% 2) Run all the tests in the repo besides tests for the `+meta/` dir
%   [status, failures] = toupee.tests.runAll('ignore', 'meta');

import toupee.misc.iif

% If we're not ignoring any folders, then run all tests.
if nargin < 1
    tests = testsuite('IncludeSubfolders', true);
    results = run(tests);
    status = iif(isempty(results) || isempty(results.Failed), true, false);
end

end

