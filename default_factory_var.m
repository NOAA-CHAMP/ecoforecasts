function result = default_factory_var(facname)
%function result = default_factory_var(facname)
%
% Return default variable name for fact-factories named FACNAME. If FACNAME
% not given or empty, return cell array of ALL default factory-variable names.  
%
% Last Saved Time-stamp: <Thu 2010-02-18 11:56:40 Eastern Standard Time Lew.Gramer>

  default_factory_var_table = { ...
      'sst_7d'			'misst_sst_7_day_maximum' ; ...
      'sst_3d'			'misst_sst_3_day_lowpass' ; ...
      'qscat_7d'		'qscat_speed_7_day_maximum' ; ...
      'global_par_7d'		'global_par_7_day_maximum' ; ...
      'par_7d'			'bic_surf_par_7_day_maximum' ; ...
                   };

  if ( ~exist('facname','var') || isempty(facname) )
    result = default_factory_var_table;

  else
    result = '';
    ix = find(strcmpi(default_factory_var_table(:,1), facname));
    if ( ~isempty(ix) )
      result = default_factory_var_table{ix,2};
    end;
  end;

return;
