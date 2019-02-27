
M1 = dlmread('result-iso-len.txt');


figure(1)
hold on
scatter(M1(:,1), log10(M1(:,2)))

title('computation time (100 nt)','fontsize',24)
xlabel('Hamming distance','fontsize',24)
ylabel('log_{10} (computation time (second))','fontsize',24)
set(gca,'FontSize',18);