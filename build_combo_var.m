function stn = build_combo_var(stn, varname)
%function stn = build_combo_var(stn, varname)
%
% A "combo variable" is a combination of two other time series, using one of
% the following operators to combine each value within each of the two ts'es:
%  combo_ops = {merge avg sum diff var prod div max min spd dir};
%
% As of 2018 Feb 25, will also build derived time series for each colummn in
% field STN.(varname).prof, if found. If so, CALLS: FILTER_GAPS_PROF (v.)
%
% Last Saved Time-stamp: <Sun 2018-02-25 16:22:55 Eastern Standard Time gramer>


  % Don't rebuild if we already have everything in place!
  if ( isfield(stn, varname) )
    return;
  end;

  combo_ops = {'merge' 'avg' 'sum' 'diff' 'var' 'prod' 'div' 'max' 'min' 'spd' 'dir'};

  for opc = combo_ops

    op = opc{:};
    findidx = regexp(varname, ...
                     ['_[0-9]+_(week|w|day|d|hour|h|minute|m|second|s)s?_[a-z]+_(' op '$|' op '_)']);
    for idx = 1:length(findidx)

      %DEBUG:      disp(['Building combo ' varname ' (' op ')']);

      begidx = findidx(idx);
      endidx = regexp(varname(begidx:end), ['_(' op '$|' op '_)']);
      endidx = begidx + endidx(1) - 1;

      % Make sure each of two constituent derived variables has been built
      var1 = varname(1:begidx-1);
      %DEBUG:      disp(['Building derived ' var1]);
      stn = build_derived_var(stn, var1);
      var2 = varname(endidx+length(op)+2:end);
      %DEBUG:      disp(['Building derived ' var2]);
      stn = build_derived_var(stn, var2);

      if ( ~isfield(stn, var1) || ~isfield(stn.(var1), 'data') || ...
           ~isnumeric(stn.(var1).data) )
        error('No valid base field #1 "%s" in stn!', var1);
      end;
      if ( ~isfield(stn, var2) || ~isfield(stn.(var2), 'data') || ...
           ~isnumeric(stn.(var2).data) )
        error('No valid base field #2 "%s" in stn!', var2);
      end;

      sep = find(varname(begidx:endidx-1) == '_', 1, 'last') + begidx - 1;
      per = varname(begidx+1:sep-1);
      derop = varname(sep+1:endidx-1);

      switch ( derop ),
       case {'asof', 'current'},
        dervar1 = var1;
        dervar2 = var2;
       otherwise,
        dervar1 = [var1 '_' per '_' derop];
        dervar2 = [var2 '_' per '_' derop];
        %DEBUG:        disp(['Building dervar ' dervar1]);
        stn = build_derived_var(stn, dervar1);
        %DEBUG:        disp(['Building dervar ' dervar2]);
        stn = build_derived_var(stn, dervar2);
      end;

      if ( ~isfield(stn, dervar1) || ~isfield(stn.(dervar1), 'data') || ...
           ~isnumeric(stn.(dervar1).data) )
        error('Unable to build combo field #1 "%s" in stn!', dervar1);
      end;
      if ( ~isfield(stn, dervar2) || ~isfield(stn.(dervar2), 'data') || ...
           ~isnumeric(stn.(dervar2).data) )
        error('Unable to build combo field #2 "%s" in stn!', dervar2);
      end;


      doProf = false;
      if ( isfield(stn.(var1),'prof') && size(stn.(var1).prof,1) == numel(stn.(var1).date) && ...
           isfield(stn.(var2),'prof') && size(stn.(var2).prof,1) == numel(stn.(var2).date) )
        doProf = true;
      end;

      %DEBUG:      disp([varname,' = ',dervar1,' ',upper(op),' ',dervar2]);

      if ( strcmp(op, 'merge') )

        stn.(varname).date = union(stn.(dervar1).date, stn.(dervar2).date);
        s1 = find(ismember(stn.(varname).date, stn.(dervar1).date));
        s2 = find(ismember(stn.(varname).date, stn.(dervar2).date));
        stn.(varname).data(s1,1) = stn.(dervar1).data;
        stn.(varname).data(s2,1) = stn.(dervar2).data;
        if ( doProf )
          stn.(varname).prof(s1,:) = stn.(dervar1).prof;
          stn.(varname).prof(s2,:) = stn.(dervar2).prof;
        end;

      else

        stn.(varname).date = intersect(stn.(dervar1).date, stn.(dervar2).date);
        s1 = find(ismember(stn.(dervar1).date, stn.(varname).date));
        s2 = find(ismember(stn.(dervar2).date, stn.(varname).date));

        v1 = stn.(dervar1).data(s1);
        v2 = stn.(dervar2).data(s2);
        if ( doProf )
          p1 = stn.(dervar1).prof(s1,:);
          p2 = stn.(dervar2).prof(s2,:);
        end;

        switch (op)
         case 'avg',
          stn.(varname).data = (v1 + v2) / 2;
          if doProf; stn.(varname).prof = (p1 + p2) / 2; end;
         case 'sum',
          stn.(varname).data = (v1 + v2);
          if doProf; stn.(varname).prof = (p1 + p2); end;
         case 'diff',
          stn.(varname).data = (v1 - v2);
          if doProf; stn.(varname).prof = (p1 - p2); end;
         case 'var',
          stn.(varname).data = abs(v1 - v2);
          if doProf; stn.(varname).prof = abs(p1 - p2); end;
         case 'prod',
          stn.(varname).data = (v1 .* v2);
          if doProf; stn.(varname).prof = (p1 .* p2); end;
         case 'div',
          stn.(varname).data = (v1 ./ v2);
          if doProf; stn.(varname).prof = (p1 ./ p2); end;
         case 'max',
          stn.(varname).data = max(v1, v2);
          if doProf; stn.(varname).prof = max(p1, p2); end;
         case 'min',
          stn.(varname).data = min(v1, v2);
          if doProf; stn.(varname).prof = min(p1, p2); end;
         case 'spd',
          stn.(varname).data = sqrt((v1.*v1) + (v2.*v2));
          if doProf; stn.(varname).prof = sqrt((p1.*p1) + (p2.*p2)); end;
         case 'dir',
          %stn.(varname).data = asind((-v1) ./ sqrt((v1.*v1) + (v2.*v2)));
          stn.(varname).data = uv_to_dir(v1, v2);
          if doProf; stn.(varname).prof = uv_to_dir(p1, p2); end;
         otherwise,
        end;

      end;

    end;

  end;

return;
