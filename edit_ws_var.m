function edit_ws_var
%function edit_ws_var
%
% Edit a variable by evaluating it in the base (global) workspace, and
% calling 'openvar' (cf. m-file help)
%
% Last Saved Time-stamp: <Wed 2008-08-13 08:23:21 Eastern Daylight Time gramer>
%

    evalin('base', 'x = { ''lo'' 0 1 ; ''av'' 1 3 ; ''hi'' 3 6 };');
    openvar('x');

    x = evalin('base', 'x');

return;
