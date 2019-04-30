1;

stanms = {'lkwf1','fwyf1','mlrf1','smkf1','sanf1','dryf1','42003'};

for ix = 1:length(stanms)
  nm = stanms{ix};

  delete([nm '-stats.csv']);

  fprintf('\n =================== %s : ALL ===================\n', nm);
  stn = anndbcs(nm,{'annual','winter','spring','summer','autumn'},[],[],[],true);
  fprintf('\n =================== %s : 93%% ===================\n', nm);
  stn = anndbcs(stn,{'annual','winter','spring','summer','autumn'},[],[],93,true);

  stn = [];
  clear stn;
end;
