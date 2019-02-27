
M1 = dlmread('result-iso-dist.txt');


figure(1)
hold on
scatter(M1(:,1), log10(M1(:,2)))

title('computation time (distance 18)','fontsize',24)
xlabel('RNA length (nt)','fontsize',24)
ylabel('log_{10} (computation time (second))','fontsize',24)
set(gca,'FontSize',18);