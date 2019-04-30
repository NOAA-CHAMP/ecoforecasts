function cod = rgb2color(clr)
%function cod = rgb2color(clr)
% Converts a HEX or RGB color code to a MATLAB RGB color code (Matlab RGB values are displayed from 0 to 1)
%   Input color formats allowed
%   HEX: '#FFFFFF' or 'FFFFFF'
%   RGB: [52 152 219]
% FROM: http://codegists.com/code/plot-color-matlab-rgb/

  inputSize = length(clr);
 
  switch ( inputSize ),
   case {6,7},  % HEX string 'FFFFFF' or '#FFFFFF'
    if  inputSize == 7
      clr = clr(2:end);
    end
    r = double(hex2dec(clr(1:2)))/255;
    g = double(hex2dec(clr(3:4)))/255;
    b = double(hex2dec(clr(5:6)))/255;
    cod = [r, g, b];

   case 3,  % RGB array [123 123 123]
    for ix = 1:3
      if ( clr(ix) < 0 || clr(ix) > 255 )
        error('RGB values must be between 0 and 255');
      end;
    end;
    r = clr(1)/255;
    g = clr(2)/255;
    b = clr(3)/255;
    cod = [r, g, b];
    
   otherwise,
    error('Argument must be either a HEX string in "#FFFFFF" or "FFFFFF" format, or an RGB array in [255 255 255] format');
  end;

return;
