function [X,ndx,dbg] = natsortrows(X,col,varargin)
% Alphanumeric / Natural-Order sort the rows of a cell array of strings (1xN char).
%
% (c) 2012 Stephen Cobeldick
%
% Alphanumeric sort of the rows of a cell array of strings: sorts by both
% character order and also the values of any numbers that occur within the
% strings. The cell of strings must be a matrix. SORTROWS input <col> is
% also supported, so NATSORTROWS is a drop-in replacement for SORTROWS.
%
%%% Example:
% natsortrows({'x2','b';'x10','a';'x1','a';'x2','a'})
%  ans = 
%    'x1'     'a'
%    'x2'     'a'
%    'x2'     'b'
%    'x10'    'a'
%
%%% Syntax:
%  Y = natsortrows(X)
%  Y = natsortrows(X,col)
%  Y = natsortrows(X,col,xpr)
%  Y = natsortrows(X,col,xpr,<options>)
% [Y,ndx] = natsortrows(X,...)
% [Y,ndx,dbg] = natsortrows(X,...)
%
% To sort filenames or filepaths correctly use NATSORTFILES (File Exchange 47434).
% To sort all of the strings in a cell array use NATSORT (File Exchange 34464).
%
% See also NATSORT NATSORTFILES SORTROWS SORT CELLSTR IREGEXP REGEXP SSCANF
%
%% File Dependency %%
%
% NATSORTROWS requires the function NATSORT (File Exchange 34464). The inputs
% <xpr> and <options> are passed directly to NATSORT: see NATSORT for case
% sensitivity, sort direction, numeric substring matching, and other options.
%
%% Examples %%
%
% A = {'B','2','X';'A','100','X';'B','10','X';'A','2','Y';'A','20','X'};
% sortrows(A) % wrong numeric order:
%  ans =
%    'A'  '100'  'X'
%    'A'    '2'  'Y'
%    'A'   '20'  'X'
%    'B'   '10'  'X'
%    'B'    '2'  'X'
% natsortrows(A) % correct numeric order:
%  ans =
%    'A'    '2'  'Y'
%    'A'   '20'  'X'
%    'A'  '100'  'X'
%    'B'    '2'  'X'
%    'B'   '10'  'X'
% natsortrows(A,'descend')
%  ans =
%    'B'   '10'  'X'
%    'B'    '2'  'X'
%    'A'  '100'  'X'
%    'A'   '20'  'X'
%    'A'    '2'  'Y'
%%% Sort ascending by the second column, descending by the third column:
% sortrows(A,[2,-3]) % wrong numeric order:
%  ans =
%    'B'   '10'  'X'
%    'A'  '100'  'X'
%    'A'    '2'  'Y'
%    'B'    '2'  'X'
%    'A'   '20'  'X'
% natsortrows(A,[2,-3]) % correct numeric order:
%  ans =
%    'A'    '2'  'Y'
%    'B'    '2'  'X'
%    'B'   '10'  'X'
%    'A'   '20'  'X'
%    'A'  '100'  'X'
%
% B = {'-12';'ABCD';'3e45';'67.8';'+9';'+Inf';'NaN'};
% sortrows(B)
%  ans =
%    '+9'
%    '+Inf'
%    '-12'
%    '3e45'
%    '67.8'
%    'ABCD'
%    'NaN'
% natsortrows(B,[],'NaN|(+|-)?(Inf|\d+\.?\d*((e|E)(+|-)?\d+)?)')
%  ans =
%    '-12'
%    '+9'
%    '67.8'
%    '3e45'
%    '+Inf'
%    'NaN'
%    'ABCD'
%
%% Input and Output Arguments %%
%
% See NATSORT for a full description of <xpr> and the <options>.
%
%%% Inputs (*=default):
%  X   = CellArrayOfCharRowVectors, size MxN. With rows to be sorted.
%  col = NumericVector, column indices to sort <X> by the corresponding
%        columns, where >0=ascending, <0=descending, exactly as per SORTROWS.
%      = CharRowVector, 'descend'/'ascend'*, sort direction selection.
%  xpr = CharRowVector, regular expression to detect numeric substrings.
%  <options> can be supplied in any order and are passed directly to NATSORT.
%            Excludes 'descend'/'ascend', which can be supplied as input <col>.
%
%%% Outputs:
%  Y   = CellArrayOfCharRowVectors, input <X> with the rows sorted as per <col>.
%  ndx = NumericVector, size Mx1. Row indices such that Y = X(ndx,:).
%  dbg = CellVectorOfCellArrays, size 1xN. Each cell contains the debug cell
%        array for one column of input <X>. To help debug <xpr>. See NATSORT.
%
% [Y,ndx,dbg] = natsortrows(X,*col,*xpr,<options>)

%% Input Wrangling %%
%
assert(iscell(X),'First input <X> must be a cell array.')
tmp = cellfun('isclass',X,'char') & cellfun('size',X,1)<2 & cellfun('ndims',X)<3;
assert(all(tmp(:)),'First input <X> must be a cell array of char row vectors (1xN char).')
assert(ismatrix(X),'First input <X> must be a matrix (size RxC).')
%
%% Select Columns to Sort %%
%
[m,n] = size(X);
dbg{n} = [];
ndx = 1:m;
drn = {'descend','ascend'};
isn = false;
%
if nargin<2 || isnumeric(col)&&isempty(col)
	vec = n:-1:1;
elseif ischar(col)&&isrow(col)&&any(strcmpi(col,drn))
	% Sort all columns descending/ascending.
	vec = n:-1:1;
	varargin{max(2,end+1)} = col;
elseif isnumeric(col)
	% Sort columns according to the provided indices.
	assert(isreal(col)&&isvector(col),'Second input <col> must be a real numeric vector.')
	assert(all(fix(col)==col)&&all(abs(col)<=n)&&all(col),...
		'Second input <col> must be a vector of column indices into the first input <X>.')
	vec = reshape(col(end:-1:1),1,[]);
	varargin{max(2,end+1)} = [];
	isn = true;
else
	error('Second input <col> must be a numeric vector of indices, or ''ascend''/''descend''.')
end
%
%% Sort Columns %%
%
for k = vec
	if isn
		varargin(end) = drn((3+sign(k))/2);
		k = abs(k); %#ok<FXSET>
	end
	if nargout<3 % faster:
		[~,ids] = natsort(X(ndx,k),varargin{:});
	else % for debugging:
		[~,ids,tmp] = natsort(X(ndx,k),varargin{:});
		[~,idd] = sort(ndx);
		dbg{k} = tmp(idd,:);
	end
	ndx = ndx(ids);
end
%
ndx = ndx(:);
X = X(ndx,:);
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%natsortrows