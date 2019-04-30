1;

% Make ALL graphs invisible - very handy for batch matlab runs with plots in
% them - and in particular, plots which we DO want to be created and saved to
% a PS or other file, but which we do NOT want to suddenly pop up all over
% our screens while we are trying to work on other things... :)
set(0, 'DefaultFigureVisible', 'off');

disp('MATLAB Graphs will *NOT* be visible by default from now on...');
