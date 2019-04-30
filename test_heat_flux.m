1;

wz = 35; az = 30;

w = 0:0.1:25;
p = 1014;
t = 25;

sf = figure; set(gcf,'units','norm','outer',[0 0 1 1]);
lf = figure; set(gcf,'units','norm','outer',[0 0 1 1]);

ls = repmat('rgbcmyk', [1 10]);

qs = 50:25:100;
qs = [50 75 90 99];
as = 20:2.5:30;

for qi = 1:length(qs)
  q = qs(qi);
  disp(q);
  for ai = 1:length(as)
    a = as(ai);

    result = hfbulktc(w,wz,a,az,q,az,p,t);
    wind_stress = result(:,4);
    sensible_heat_flux = result(:,1);
    latent_heat_flux = result(:,2);

    figure(sf); subplot(1,length(qs),qi); hold on; plot(w,sensible_heat_flux,ls(ai));
    figure(lf); subplot(1,length(qs),qi); hold on; plot(w,latent_heat_flux,ls(ai));
    lgnd{ai} = ['T_a_i_r=' num2str(a)];
  end;
  figure(sf); subplot(1,length(qs),qi); legend(lgnd,'Location','Best'); title(['RH=' num2str(q)]);
  figure(lf); subplot(1,length(qs),qi); legend(lgnd,'Location','Best'); title(['RH=' num2str(q)]);
end;

figure(sf); suptitle('COARE 2.0: H_S_e_n_s_i_b_l_e');
figure(lf); suptitle('COARE 2.0: H_L_a_t_e_n_t');
