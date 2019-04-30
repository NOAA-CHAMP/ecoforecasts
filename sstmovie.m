1;

for j=1:n
  plot_command
  M(j) = getframe;
end
movie(M);
movie2avi(M,'sstmovie.avi');
