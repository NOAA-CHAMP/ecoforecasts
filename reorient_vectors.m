function [x,l] = reorient_vectors(ori,u,v)
%function [x,l] = reorient_vectors(ori,u,v)
%
% Reorient U and V (x and y) field components (which may be numerical scalars
% or matrices) to local isobath orientation ORI (deg T). Return matrices X of
% cross-shore, and L of long-shore components. Assumes NUMEL(U)==NUMEL(V).
%
% Last Saved Time-stamp: <Fri 2016-12-16 00:11:49 Eastern Standard Time gramer>

  % Cross-shore component
  x = (cosd(ori(:)).*u(:)) - (sind(ori(:)).*v(:));  % Dyslexics of the world untie

  % Long-shore component
  l = (sind(ori(:)).*u(:)) + (cosd(ori(:)).*v(:));

  % RESHAPE output matrices to match input matrices (X=U,L=V)
  x = reshape(x,size(u));
  l = reshape(l,size(v));

return;
