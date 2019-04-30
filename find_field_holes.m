function edgeix = find_field_holes(d,varargin)
%function edgeix = find_field_holes(d[,holeFunc][,bdrySize,useEightConnected,useEdges])
%
% Return edge indices of each "hole" (group of contiguous NaNs) in field D.
% Specify FOURNOTEIGHT true if corners should be ignored. Specify USEEDGES if
% the edges of the field (D(:,[1,end]) and D([1,end],:) should be considered.
%
% Last Saved Time-stamp: <Sat 2016-05-14 13:56:29 Eastern Daylight Time gramer>

  if ( ~isnumeric(d) || ndims(d) ~= 2 )
    error('First arg should be a 2D numeric matrix');
  end;

  args = varargin;

  holeFunc = @isnan;
  if ( numel(args) > 1 && isa(args{1},'function_handle') )
    holeFunc = args{1};
    args(1) = [];
  end;
  bdrySize = 1;
  if ( numel(args) >= 1 && ~isempty(args{1}) )
    bdrySize = args{1};
    args(1) = [];
  end;
  useEightConnected = false;
  if ( numel(args) >= 1 && ~isempty(args{1}) )
    useEightConnected = args{1};
    args(1) = [];
  end;
  useEdges = false;
  if ( numel(args) >= 1 && ~isempty(args{1}) )
    useEdges = args{1};
    args(1) = [];
  end;

  % Find all indices in "holes" (NaNs)
  ix = find(holeFunc(d));
  [yix,xix] = ind2sub(size(d),ix);


  edgeix = [];

  % Remove field edges
  if ( useEdges )
    edgeix = ix(yix<=bdrySize | xix<=bdrySize | yix>size(d,1)-bdrySize | xix>size(d,2)-bdrySize);
  end;

  ix(yix<=bdrySize | xix<=bdrySize | yix>size(d,1)-bdrySize | xix>size(d,2)-bdrySize) = [];
  [yix,xix] = ind2sub(size(d),ix);

  if ( ~useEightConnected )
    for dx=-bdrySize:bdrySize;
      nix = sub2ind(size(d),yix,xix+dx);
      nbix = find(~isnan(d(nix)));
      edgeix(end+1:end+numel(nbix)) = nix(nbix);
    end;
    for dy=-bdrySize:bdrySize;
      nix = sub2ind(size(d),yix+dy,xix);
      nbix = find(~isnan(d(nix)));
      edgeix(end+1:end+numel(nbix)) = nix(nbix);
    end;

  else
    for dx=-bdrySize:bdrySize;
      for dy=-bdrySize:bdrySize;
        if (dx || dy)
          nix = sub2ind(size(d),yix+dy,xix+dx);
          nbix = find(~isnan(d(nix)));
          edgeix(end+1:end+numel(nbix)) = nix(nbix);
        end;
      end;
    end;
  end;

return;
