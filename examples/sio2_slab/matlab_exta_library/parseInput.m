function inst_out = parseInput(paramMatch, varargin)
% *************************************************************************
% Parse input arguments from varagin based on instructions on paramMatch
% The instructions are stored in the structure "inst_out"
% Input:
%   paramMatch: cell array instruction name's and default values
%   varargin: cell array with input instructions
% Output:
%   inst_out: structure with the set instructions
%
% F. Ramirez 2020
% *************************************************************************

p = inputParser();
p.KeepUnmatched = false;
p.CaseSensitive = false;
p.StructExpand  = false;
p.PartialMatching = true;

for i = 1:2:length(paramMatch)
    addParameter( p, paramMatch{i}  , paramMatch{i+1} )
end
parse( p, varargin{:} )

inst_out = p.Results;
end