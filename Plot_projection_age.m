
% load data in the form of calendar age (year BP), 2sd, 14C age (year BP), 2sd


%Burdwood Bank deep-sea coral data at water depth of 1879m as an example 
data = [797	9	2199	52
2108	93	3316	54
2139	40	3342	54
2530	31	3623	54
2684	40	3838	51
2715	37	3814	52
5390	81	6063	55
5556	43	6216	54
6063	48	6621	59
6843	82	7427	56];

% Calculate the projection ages and error ellipses using the Projection_Age
% function
[meanage,meanProjVentAge,meanProjDD14Ccorr,ProjVentAge_ellipse, ...
    ProjDD14Ccorr_ellipse] = Projection_Age(data);

%plot projection age results
subplot(111);hold on
for ii=1:length(meanage)
    plot(ProjVentAge_ellipse(:,ii*2-1),ProjVentAge_ellipse(:,ii*2),'k-');
end
plot(meanage,meanProjVentAge,'ko','markersize',7);
xlabel('Age (year BP)');
ylabel('Projection age (year)');