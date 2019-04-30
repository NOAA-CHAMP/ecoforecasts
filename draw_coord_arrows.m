function draw_coord_arrows(varargin)
%function draw_coord_arrows(varargin)
%
% function out = draw_arrow(basept,headpt,headsize) 
% %by Ryan Molecke


%   basept = varargin{1};
%   len = varargin{2};
%   ang = varargin{3};

%   headpt = basept + 

%   v1 = len/2.5;


   basept = varargin{1};
   headpt = varargin{2};
   if ( nargin > 2 )
     headsize = varargin{3};
   else
     headsize = 0.5;
   end;

   v1 = headsize*(basept-headpt)/2.5;

   theta = 22.5*pi/180; 
   theta1 = -1*22.5*pi/180; 
   rotMatrix = [cos(theta) -sin(theta) ; sin(theta) cos(theta)]; 
   rotMatrix1 = [cos(theta1) -sin(theta1) ; sin(theta1) cos(theta1)]; 
    
   v2 = v1*rotMatrix; 
   v3 = v1*rotMatrix1; 
   x1 = headpt; 
   x2 = x1 + v2; 
   x3 = x1 + v3; 

    % below line fills the arrowhead (black) 
   fill([x1(1) x2(1) x3(1)],[x1(2) x2(2) x3(2)],[0 0 0]); 
   % below line draws line (black) 
   plot([basept(1) headpt(1)],[basept(2) headpt(2)],'linewidth',2,'color',[0 0 0]);


return;
