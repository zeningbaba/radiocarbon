% load data in the form of U-Th age, 2sd, 14C age, 2sd
data = [797	9	2199	52
2108	93	3316	54
2139	40	3342	54
2530	31	3623	54
2684	40	3838	51
2715	37	3814	52
5390	81	6063	55
5556	43	6216	54
6063	48	6621	59
6843	82	7427	56];%Burdwood Bank 1879m

% Calculate the delta14C values and error ellipses using the Magiccarbon_intensity
% function
[meanage,smplnumber,meandelta14C,meandeltadelta14C,deltaellipse, ...
    deltadeltaellipse] = Delta14C(data,0);


%plot delta14C results
subplot(111);hold on
for ii=1:length(meanage)
    plot(deltaellipse(:,ii*2-1),deltaellipse(:,ii*2),'k-');
end
plot(meanage,meandelta14C,'ko','markersize',7);
xlabel('Age (year BP)');
ylabel('Delta14C');

