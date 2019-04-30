function nil = empty_ts
%function nil = empty_ts
% Return an empty Time Series struct, i.e., with empty .date and .data fields.
  nil=struct('date',[],'data',[]);
return;
