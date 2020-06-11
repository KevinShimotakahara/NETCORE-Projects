figure(2)
h = gca;
subplot(4, 1, 1);
hold on
grid on
%xlabel('Time (s)','FontSize',16)
ylabel('Power (W)','FontSize',16)
%ylim([0 9000])
plot(ScopeSPS4.time,ScopeSPS4.signals(3).values+ScopeSPS4.signals(4).values+ScopeSPS4.signals(5).values,'m')
plot(ScopeSPS.time,ScopeSPS.signals(3).values+ScopeSPS.signals(4).values+ScopeSPS.signals(5).values,'r')
plot(ScopeSPS2.time,ScopeSPS2.signals(3).values+ScopeSPS2.signals(4).values+ScopeSPS2.signals(5).values,'b')
plot(ScopeSPS3.time,ScopeSPS3.signals(3).values+ScopeSPS3.signals(4).values+ScopeSPS3.signals(5).values,'k')

title('Demand Response Household Aggregate Demand')
legend('No DR','LTE R12 D2D TM-2', 'HRQ','Ideal Data Transfer')
h = gca;
set(h,'FontSize',14,'FontWeight','Bold')

subplot(4,1,2)
hold on
grid on
%xlabel('Time (s)','FontSize',16)
ylabel('Power (W)','FontSize',16)
yticks(-10000:5000:10000)
plot(ScopeSPS4.time,-ScopeSPS4.signals(2).values,'m')
plot(ScopeSPS.time,-ScopeSPS.signals(2).values,'r')
plot(ScopeSPS2.time,-ScopeSPS2.signals(2).values,'b')
plot(ScopeSPS3.time,-ScopeSPS3.signals(2).values,'k')

title('Transformer Power Flow')
leg = legend('No DR','LTE R12 D2D TM-2', 'HRQ','Ideal Data Transfer');
leg.Location = 'southeast';
h = gca;
set(h,'FontSize',14,'FontWeight','Bold')
h.YRuler.Exponent = 0;

subplot(4,1,3)
hold on
grid on
%xlabel('Time (s)','FontSize',16)
ylabel('Power (W)','FontSize',16)
%ylim([-6000 9000])
plot(ScopeSPS.time,ScopeSPS.signals(1).values,'r')
title('PV Power Generation')
h = gca;
set(h,'FontSize',14,'FontWeight','Bold')

subplot(4,1,4)
hold on
grid on
title('Electricity Price Curve')
xlabel('Time (s)','FontSize',16)
ylabel('Price (Arbitrary Units)','FontSize',16)
plot(ScopeSPS.time,ScopeSPS.signals(6).values,'b')
legend('Price')
h = gca;
set(h,'FontSize',14,'FontWeight','Bold')