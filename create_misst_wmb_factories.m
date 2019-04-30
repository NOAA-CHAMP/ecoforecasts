function stns = create_misst_wmb_factories(stns)

  [stns.factories] = deal([]);

  for ix = 1:length(stns)
    PRCTILES = [];

    stns(ix).factories.seatemp_7d = ...
        create_prctile_factory(stns(ix), 'misst_sst_7_day_maximum', PRCTILES);

    stns(ix).factories.wind_7d = ...
        create_prctile_factory(stns(ix), 'qscat_speed_7_day_maximum', PRCTILES);

    stns(ix).factories.par_7d = ...
        create_prctile_factory(stns(ix), 'nesdis_satpar_7_day_maximum', PRCTILES);
  end;

return;
