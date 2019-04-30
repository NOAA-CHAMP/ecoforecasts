function U=underline(S,fns) 
if nargin<2 
    fns = 10; 
end 
p=''; 
for n = 1:length(S) 
    p=[p '\_\_']; 
end 
U =['_{\fontsize{' num2str(fns) '}^{' S '}_{^{' p '}}}']; 
return;
