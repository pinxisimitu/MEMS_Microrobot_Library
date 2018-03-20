function out = etch_points(h_etch,slice)
% This function helps generate the points for the etch hole generation
% script. 
% Needs h_etch parameters of etch hole radius and undercut

% Returns a vector (or single point) that is the offest from the inital
% point (p0) at which the etch hole array should start. This represents the
% center of the etch hole.

%slice is the width or length that needs to be populated with etch holes

% Determine where the etch hole offset are 

if slice < (2*(h_etch.undercut + h_etch.undercut/sqrt(2)) + 2*h_etch.r)
    % Only a single row/column
    out = slice/2;
elseif slice < (2*(h_etch.undercut + h_etch.undercut/sqrt(2)) + 4*h_etch.r + 2*h_etch.undercut/sqrt(2))
    % Two rows/columns
    gap = (slice - 2*(2*h_etch.r))/(3);
    out = (gap + 2*h_etch.r)*[1:2]-h_etch.r;
elseif slice < (2*(h_etch.undercut + h_etch.undercut/sqrt(2)) + 6*h_etch.r + 4*h_etch.undercut/sqrt(2))
    % Three rows/columns
    gap = (slice - 3*(2*h_etch.r))/(4);
    out = (gap + 2*h_etch.r)*[1:3]-h_etch.r;
else
    % There are more than three rows/columns
    n = round((slice - 2*(h_etch.undercut + h_etch.undercut/sqrt(2)) + 2*h_etch.undercut/sqrt(2))/(2*h_etch.r+(2*h_etch.undercut/sqrt(2))));
    gap = (slice - n*(2*h_etch.r))/(n+1);
    out = (gap + 2*h_etch.r)*[1:n]-h_etch.r;
end
