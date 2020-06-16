function [status, results] = runToupeeTests(varargin)
% Runs all tests in the 'toupee' repo
% Gathers and runs all tests in the recursive contents of this directory.
% All input arguments should be name-value pairs.
% 
% Inputs:
% -------
% 'ignore' (optional) : cell
%   The folder directories of tests to *not* include.      
% 
% Outputs:
% --------
% status : logical
%   0 if all the tests that were run passed, 1 otherwise.
% results : TestResult
%   Full details on all the tests that were run.
%
% Examples:
% ---------
% 1) Run all the tests in the repo:
%   [status, failures] = runToupeeTests();
% 2) Run all the tests in the repo besides tests for the `+meta/` dir
%   [status, failures] = runToupeeTests('ignore', 'meta');

import toupee.misc.iif

% Navigate to this function's directory for running tests, and return to
% the current MATLAB working directory when this function ends
cwd = pwd();
[twd, ~, ~] = fileparts(mfilename('fullpath'));
cd(twd);
c = onCleanup(@() cd(cwd));

% If we're not ignoring any folders, then run all tests.
if nargin < 1
    tests = testsuite('IncludeSubfolders', true);
    results = run(tests);
    status = iif(isempty(results) || isempty(results.Failed), true, false);
end

end

