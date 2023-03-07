
clear
clf

voltage_range = (-120 :.1: 40);

tauH_out = zeros(length(voltage_range),1);
tauM_out = zeros(length(voltage_range),1);
Hinfintiy_out = zeros(length(voltage_range),1);
Minfinity_out = zeros(length(voltage_range),1);

for i= 1:length(voltage_range)
    clear tauH tauM Minfinity Hinfinity
    setGlobal2(voltage_range(i))
    hh_core()
    global tauH tauM Minfinity Hinfinity
    tauH_out(i)= tauH;
    tauM_out(i) = tauM;
    Hinfintiy_out(i)= Hinfinity;
    Minfinity_out(i) = Minfinity;
end



%% Ploting 

figure(1);
subplot(2,1,1)
hold on
yyaxis left
plot(voltage_range, tauM_out)
yyaxis right
plot(voltage_range,tauH_out)
xlabel('voltage')
ylabel('tau')
legend('tauM','tauH');
hold off

subplot(2,1,2)
hold on
yyaxis left
plot(voltage_range,Minfinity_out)
yyaxis right
plot(voltage_range, Hinfintiy_out)
xlabel('voltage')
ylabel('Infinities')
legend('Minfinity','Hinfinity');
    
function setGlobal2(voltage)
global V_insert
V_insert = voltage;
end

