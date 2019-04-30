1;

x = 0 + 0:(1/24):31;

% tide_comps = [ ...
%              ];

aK1 = 8.3/1e2;	pK1 = (23.934/24)/(2*pi);	phK1 = (-172/180)*pi;
aO1 = 8.2/1e2;	pO1 = (25.819/24)/(2*pi);	phO1 = (-39/180)*pi;

aM2 = 17.8/1e2;	pM2 = (12.42/24)/(2*pi);	phM2 = (-15/180)*pi;
aS2 = 5.1/1e2;	pS2 = (12.00/24)/(2*pi);	phS2 = (149/180)*pi;

amp = (aK1.*sin((x/pK1)+phK1)) + ...
      (aO1.*sin((x/pO1)+phO1)) + ...
      (aM2.*sin((x/pM2)+phM2)) + ...
      (aS2.*sin((x/pS2)+phS2));

figure;
plot(x, amp);
