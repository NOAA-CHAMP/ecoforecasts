function varargout = colorbar_tight(varargin)
%function varargout = colorbar_tight(varargin)
%
% CALLS: varargout{:} = colorbar(varargin{:}); tighten_axes;

  varargout{:} = colorbar(varargin{:});
  tighten_axes;

return;
