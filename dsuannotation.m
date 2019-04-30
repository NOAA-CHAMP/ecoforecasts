function anh = dsuannotation(antyp,varargin)
%function anh = dsuannotation(antyp,[X,Y | POS][,ANNOTATION_ARGS])
%
% Add ANNOTATION (v.) of type ANTYP to current FIGURE, accepting coordinates
% (X and Y, or POS 4-vector) in data space units (e.g., longitude, latitude)
% rather than "normalized" (0-1) figure units. X,Y may be Nx2 matrices of
% data space units or they may be Nx3 (DOUBLE) matrices with rows containing
% an AXES handle, X data coordinate, and Y data coordinate. This latter form
% allows annotations to cross between AXES. Similarly, the rows of Nx5 POS
% may be [X,Y,W,H] or [AX,X,Y,W,H]. If no AXES are specified, uses GCA.
%
% CALLS: DS2NFU (convert data space coordinates into Normalized Figure Units)
%
% Last Saved Time-stamp: <Wed 2018-10-17 13:29:33 EDT lew.gramer>
  
  anh = [];
  if ( ~ischar(antyp) )
    error('First arg must be an ANNOTATION (v.) TYPE name string');
  end;
  args = varargin;
  
  arg1 = args{1};
  args(1) = [];
  nannots = size(arg1,1);
  ncols = size(arg1,2);
  
  % Some of the calculations below depend on AXES already being where
  % they are intended to ultimately be. DRAWNOW should ensure that...
  drawnow;
  
  if ( numel(args) > 0 && isnumeric(args{1}) )
    % X,Y
    arg2 = args{1};
    args(1) = [];
    if ( ndims(arg1)~=ndims(arg2) || any(size(arg1)~=size(arg2)) )
      error('Args X and Y must have the same shape!');
    end;
    switch ( ncols )
     case 2, % [X1,X2]
      xax = repmat(gca,[nannots,1]);
      yax = xax;
      xen = arg1;
      yen = arg2;
     case 4, % [AX1,X1,AX2,X2] - AX1,AX2 (numerical) AXES handles
      xax = arg1(:,[1,3]);
      yax = arg2(:,[1,3]);
      xen = arg1(:,[2,4]);
      yen = arg2(:,[2,4]);
     otherwise,
      error('If X,Y are given, each must be either an Nx2 or Nx4 matrix');
    end;
  else
    % POS
    switch ( ncols )
     case 4, % Nx4 position matrices POS
      ax = repmat(gca,[nannots,1]);
      pos = arg1;
     case 5, % [AX,POS] position matrices POS including (numerical) AXES handles
      ax = arg1(:,1);
      pos = arg1(:,2:end);
     otherwise,
      error('If POS is given, each must be either an Nx4 or Nx5 matrix');
    end;
  end;
  
  for ix=1:nannots
    if ( exist('xen','var') )
      [xen1,yen1] = ds2nfu(xax(ix,1),xen(ix,1),yen(ix,1));
      [xen2,yen2] = ds2nfu(xax(ix,2),xen(ix,2),yen(ix,2));
      anh(ix) = annotation(antyp,[xen1,xen2],[yen1,yen2],args{:});
    elseif ( exist('pos','var') )
      posn = ds2nfu(ax(ix),pos(ix,:));
      anh(ix) = annotation(antyp,posn,args{:});
    else
      error('Invalid second argument??');
    end;
  end;
  
return;
