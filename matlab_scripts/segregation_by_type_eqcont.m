function [S_all, S_same, S_other, W_same, B_all, B_same, B_other] = segregation_by_type_eqcont( M, Ci, Ti )
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% segregation_by_type_eqcont.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION:
%    Calculate the segregation of groups of systems (i.e., system type) 
%    in a way where ontribution of each 'system' is the same regardless of 
%    its size. This is the version used in Chan et al. 2014 PNAS.
%
% USAGE:
%    [S_all, S_same, S_other, W_same, B_all, B_same, B_other] = segregation_by_type_eqcont( M, Ci, Ti );
%
% Inputs: M   :   Weighted symmetrical matrix (n x n matrix)
%                   e.g., fisher-z-transformed correlation matrix.
%
%         Ci  :   Community affiliation vector
%                   e.g., Power labels for Chan et al 441 nodes:
%
%         Ti  :   System-type affiliation vector containing T unique
%                 'types'; a nodes with a system-type label of '0' is
%                 ignored when calculating segregation from communities of
%                 other types (S_other).
%                   e.g., '1' for Sensory-motor, '2' for Association
%
% Outputs: S_all:   for each system-type, the average segregation from all
%                   other communities that have an assigned system-type
%                       e.g., average segregation of association systems to
%                       other association systems and sensory-motor system)
%
%          S_same:  for each system-type, the average segregation from
%                   communities of the same type
%                       e.g., average association-to-association system
%                       segregation.
%
%          S_other: for each system-type, the average segregation from
%                   communities of all other types
%                       e.g., average association-to-sensory-motor system
%                       segregation.
%
%          W_same:  Average of the mean correlation between nodes within the same
%                       community of one system type
%
%          B_all:     Average of edges between nodes from different communities
%                       in all system types 
%
%          B_same:    Average of edges between nodes from different communities
%                       in the same system type
%
%          B_other:   Average of edges between nodes from different communities
%                    	  in the other system type
%
%    Note: Values in output vectors are ordered in ascending order of
%          system-type relative to the input vector Ti.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Reference: Chan et al. (2014) PNAS E4997
%   2015-2019
%   Neil Savalia, UTD
%   Micaela Chan, UTD
%
%   Modification History:
%   Dec 2015: original (NKS, MYC)
%   Oct 2018: took out warning for Inf diagonal (diagonal not used in
%             calculation; edited comments (MYC)
%   Sep 2019: Add W_same, B_all/same/other outputs to match R function
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Check inputs
 if size(Ci, 1) ~= size(M, 1) || size(Ci, 1) ~= size(M, 2) > 0
     error('Length of community label index does not match with matrix dimension.');
 elseif size(Ti, 1) ~= size(M, 1) || size(Ti, 1) ~= size(M, 2) > 0
     error('Length of system-type label index does not match with matrix dimension.');
 elseif sum(sum(sum(isnan(M)))) > 0
     disp('Warning: NaN values detected in input data matrix, may result in some NaN segregation values.');
 end

 % Index of node affiliations (community & system-type)
 nCi = unique(Ci);
 nTi = unique(Ti(Ti > 0)); % do not include system-types labeled '0'

 S_all = zeros(size(nTi));
 S_same = zeros(size(nTi));
 S_other = zeros(size(nTi));
 
 W_same = zeros(size(nTi));
 B_all = zeros(size(nTi));
 B_same = zeros(size(nTi));
 B_other = zeros(size(nTi));

 % Calculate segregation
 for x = 1:size(nTi, 1) % for each system-type

     % get communities affiliations based on current system-type
     t = unique(Ci(Ti == nTi(x))); % list of communities in current type
     ot = unique(Ci(ismember(Ti, nTi(~ismember(nTi, nTi(x)))))); % list of communities of other type(s)

     % reset within and between values
     tM = [];
     w = [];
     b.all = [];
     b.same = [];
     b.other = [];

     % for each community in current type, get pieces to calculate seg
     for y = 1:size(t, 1)

         % node indices given the current community
         sys = ismember(Ci, nCi(nCi == t(y))); % nodes in current community
         same = ismember(Ci, t(t ~= (t(y)))); % nodes in communities of same 'type' as the current community
         other = ismember(Ci, ot); % nodes in communities of other 'type(s)' relative to current community

         % within and between values
         tM = M(sys, sys);
         w(y) = mean(tM(logical(triu(ones(size(tM)), 1))));
         b.all(y) = mean(reshape(M(sys, logical(same + other)),[],1));
         b.same(y) = mean(reshape(M(sys, same),[],1));
         b.other(y) = mean(reshape(M(sys, other),[],1));

     end

     % calculate segregation
     % communities with a single node will return NaN segregation values
     
     W_same(x) = mean(w, 2);
     B_all(x) = mean(b.all, 2);
     B_same(x) = mean(b.same, 2);
     B_other(x) = mean(b.other, 2);
     
     S_all(x) = (W_same(x) - B_all(x)) ./ W_same(x);        % subject x system-type average seg to all other communities
     S_same(x) = (W_same(x) - B_same(x)) ./ W_same(x);      % subject x system-type average seg to communities of the same type
     S_other(x) = (W_same(x) - B_other(x)) ./ W_same(x);    % subject x system-type average seg to communities of other type(s)

 end

end
