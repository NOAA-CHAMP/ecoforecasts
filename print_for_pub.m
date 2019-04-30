function print_for_pub(basenm,varargin)
%function print_for_pub([fh,]basenm,[...])
% Print two versions of FIGURE FH (DEFAULT: GCF), one for preview (in PNG
% format) and one for publication (in high-resolution TIFF format). BASENM
% must be a pathname (fully qualified or otherwise) excluding file extention,
% to which both '.png' and '.tif' are appended in calls to PRINT (v.)
  
  print('-dpng',[basenm,'.png'],varargin{:});
  print('-dtiff',[basenm,'.tif'],varargin{:});

return;
