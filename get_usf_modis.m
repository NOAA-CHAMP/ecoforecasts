function [fld,lat,lon] = get_usf_modis(dtrng,hrs,rgn,typ,doPNG)
%function [fld,lat,lon] = get_usf_modis(dtrng,hrs,rgn,typ,doPNG)
%
% Download if necessary all USF Optics Lab HDF data files during the period
% DTRNG and hours HRS in region RGN, then extract all data of TYP from them.
%

  opticspath = get_ecoforecasts_path('data/usf/optics');
  % Do we want to limit ourselves to certain hours of each day?
  if ( ~exist('hrs','var') || isempty(hrs) )
    hrs = 0:24;
  end;
  if ( ~exist('rgn','var') || isempty(rgn) )
    rgn = 'SE_FL';
  end;
  if ( ~exist('typ','var') || isempty(typ) )
    typ = 'SST';
  end;
  typ = upper(typ);

  if ( ~exist('doPNG','var') || isempty(doPNG) )
    doPNG = false;
  end;

  ds = unique(floor(dtrng));

  for d = ds(:)'
    if ( d > floor(now) )
      warning('I cannot predict the future!');
      %%%% EARLY LOOP BREAK
      break;
    end;

    %http://optics.marine.usf.edu/subscription/modis/SE_FL/2012/daily/222/A20122220630.QKM.SE_FL.PASS.L3D.SST.400.png
    %http://optics.marine.usf.edu/subscription/modis/SE_FL/2012/daily/222/A20122220630.QKM.SE_FL.PASS.L3D.SST.png
    %http://optics.marine.usf.edu/subscription/modis/SE_FL/2012/daily/222/A20122220630.QKM.SE_FL.PASS.L3D.hdf
    %http://optics.marine.usf.edu/subscription/modis/SE_FL/2012/daily/222/A20122221910.QKM.SE_FL.PASS.L3D.hdf
    yr = get_year(d);
    jd = get_jday(d);

    % Directories from today or recent days may only be partial - delete them when we're done
    delDir = false;
    if ( d >= floor(now)-3 )
      delDir = true;
    end;

    baseurl = sprintf('http://optics.marine.usf.edu/subscription/modis/%s/%04d/daily/%03d/',rgn,yr,jd);

    dirfname = fullfile(opticspath,sprintf('modis.%s.%04d.daily.%03d.html',rgn,yr,jd));
    if ( ~exist(dirfname,'file') )
      urlwrite(baseurl,dirfname);
    end;
    if ( ~exist(dirfname,'file') )
      warning('Unable to download directory %s',baseurl);
      %%%% EARLY LOOP CONTINUE
      continue;
    end;

    fid = fopen(dirfname,'r');
    % html = fread(fid,Inf,'uint8=>char')';
    % fclose(fid);
    if ( fid < 3 )
      warning('Unable to read directory %s',baseurl);
      %%%% EARLY LOOP CONTINUE
      continue;
    else
      html = fread(fid,Inf,'uint8=>char')';
      fclose(fid);
    end;

    %<tr><td valign="top"><img src="/icons/unknown.gif" alt="[   ]"></td><td><a href="A20132730715.QKM.SE_FL.PASS.L3D.hdf">A20132730715.QKM.SE_FL.PASS.L3D.hdf</a></td><td align="right">18-Oct-2013 20:13  </td><td align="right">1.3M</td></tr>
    % *VS.*
    %<tr><td valign="top"><img src="/icons/image2.gif" alt="[IMG]"></td><td><a href="A20132730715.QKM.SE_FL.PASS.L3D.SST.200.png">A20132730715.QKM.SE_FL.PASS.L3D.SST.200.png</a></td><td align="right">18-Oct-2013 20:15  </td><td align="right"> 31K</td></tr>
    %<tr><td valign="top"><img src="/icons/unknown.gif" alt="[   ]"></td><td><a href="A20132730715.SE_FL.MOD35.hdf">A20132730715.SE_FL.MOD35.hdf</a></td><td align="right">18-Oct-2013 20:12  </td><td align="right"> 18K</td></tr>
    %<tr><td valign="top"><img src="/icons/image2.gif" alt="[IMG]"></td><td><a href="A20132731815.QKM.SE_FL.PASS.L3D_RRC.CI.200.png">A20132731815.QKM.SE_FL.PASS.L3D_RRC.CI.200.png</a></td><td align="right">18-Oct-2013 20:14  </td><td align="right"> 32K</td></tr>

    if ( doPNG )
      res = regexp(html, ['href="(?<fname>[AT][12][089][^"]*',typ,'[.]png)"'], 'names');
    else
      res = regexp(html, ['href="(?<fname>[AT][12][089][^"]*[.]L3D[.]hdf)"'], 'names');
    end;
    for fnix = 1:numel(res)
      fname = res(fnix).fname;
      urlfname = [baseurl,fname];
      fullfname = fullfile(opticspath,fname);
      dtres = regexp(fname, '^(?<sat>[AT])(?<yr>[12][0-9][0-9][0-9])(?<jd>[0-9][0-9][0-9])(?<hr>[0-9][0-9])(?<mn>[0-9][0-9])[.]QKM', 'names');
      hr = str2num(dtres.hr);

      if ( ismember(hr,hrs) )
        if ( ~exist(fullfname,'file') )
          %DEBUG:
          disp(['Downloading ',urlfname]);
          urlwrite(urlfname,fullfname);
        end;
        if ( ~exist(fullfname,'file') )
          warning('Unable to download file %s',urlfname);
        else
          if ( doPNG )
            fld = read_usf_modis_png(fullfname,lower(typ));
          else
            %[flds,lon,lat] = read_usf_modis_hdf(fullfname,lower(typ));
            fld = read_usf_modis_hdf(fullfname,lower(typ));
          end;
        end;
      end;

    end;

    if ( delDir )
      delete(dirfname);
    end;
  end;

return;
