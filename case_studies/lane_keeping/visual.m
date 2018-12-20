close all;
fig3 = figure;
abs_model.plot(fig3,0,'k');
abs_model.plot(fig3,1:abs_model.n,'k');
abs_model.plot_box(fig3,Safe_region,'m');
abs_model.plot(fig3,1,'c');
xlabel('v_x');
ylabel('y')
