function res = read_json(fname)

  str = importdata(fname);
  res = jsondecode(str);

return;
