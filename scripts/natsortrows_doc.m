%% NATSORTROWS Examples
% The function <https://www.mathworks.com/matlabcentral/fileexchange/47433
% |NATSORTROWS|> sorts the rows of a cell array of strings (1xN char), taking
% into account any number values within the strings. This is known as a
% _natural order sort_ or an _alphanumeric sort_. Note that MATLAB's inbuilt
% <http://www.mathworks.com/help/matlab/ref/sortrows.html |SORTROWS|>
% function sorts only by character order.
%
% For sorting filenames or filepaths use
% <https://www.mathworks.com/matlabcentral/fileexchange/47434 |NATSORTFILES|>.
%
% For sorting a cell array of strings use
% <https://www.mathworks.com/matlabcentral/fileexchange/34464 |NATSORT|>.
%
%% Basic Usage:
% By default |NATSORTROWS| interprets consecutive digits as being part of
% a single integer, each number is considered to be as wide as one letter:
A = {'B1','X2';'A1','X100';'B1','X10';'A2','X2';'A1','X20';'A10','X1';'A2','X0'};
sortrows(A)
natsortrows(A)
%% Output 2: Sort Index
% The second output argument is a numeric array of the sort indices |ndx|,
% such that |Y = X(ndx,:)| where  for |Y = natsortrows(X)|:
[~,ndx] = natsortrows(A)
%% Output 3: Debugging Array
% The third output is a cell vector of cell arrays, where each cell array
% contains individual characters and numbers (after converting to numeric).
% This is useful for confirming that the numbers are being correctly
% identified by the regular expression. The cells of the cell vector
% correspond to the columns of the original input cell matrix.
[~,~,dbg] = natsortrows(A);
dbg{:}
%% Sort Direction: |'ascend'| or |'descend'| Argument
% The second input argument may be either |'ascend'| or |'descend'|, and
% all columns will be sorted accordingly:
natsortrows(A,'ascend')
natsortrows(A,'descend')
%% Sort Direction: |SORTROWS| Column Argument
% The second input argument may be a numeric vector of column indices,
% exactly as per MATLAB's |SORTROWS|, where a positive integer will sort
% the corresponding column in ascending order, and a negative integer will
% sort the corresponding column in descending order. In this example the
% second column is sorted ascending, and the third descending:
sortrows(A,[-2,1]) % wrong numeric order:
natsortrows(A,[-2,1]) % correct numeric order:
%% Regular Expression: Decimal Numbers, E-notation, +/- Sign.
% |NATSORTROWS| is a wrapper for |NATSORT|, which means all of |NATSORT|'s
% options are also supported. In particular the number recognition can be
% customized to detect numbers with decimal digits, E-notation, a +/- sign,
% or other specific features. This detection is defined by providing an
% appropriate regular expression: see |NATSORT| for details and examples.
G = {'v10.2','b'; 'v2.5','b'; 'v2.40','a'; 'v1.9','b'};
natsortrows(G) % integers, e.g. version numbers
natsortrows(G,[],'\d+\.?\d*') % decimal values
%% Regular Expression: Interactive Regular Expression Tool
% Regular expressions are powerful and compact, but getting them right is
% not always easy. One assistance is to download my interactive tool
% <https://www.mathworks.com/matlabcentral/fileexchange/48930 |IREGEXP|>,
% which lets you quickly try different regular expressions and see all of
% <https://www.mathworks.com/help/matlab/ref/regexp.html |REGEXP|>'s
% outputs displayed and updated as you type.