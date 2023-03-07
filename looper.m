%function looper(duration0, time0,Pulsing_Current_Min)
clear
clf
clc

Injection_range = (0.040 : 0.001 : 0.080)';
duration0 = 500;
time0 = 0;

% generate additional brief; pulse of injected current
%time0=20;
%duration0=3;
%Iinj0=0.0; % in mA/mm^2


frequency  = zeros(length(Injection_range), 1);

for i = 1:length(Injection_range)
    clc
    setGlobalx(Injection_range(i), duration0, time0)
    period_output = Main_hh_core();
    if ~isempty(period_output)
        %if Injection_range(i) >= Pulsing_Current_Min
        frequency(i) = 1/period_output(length(period_output));
        %print(frequency(i))
        %end
    end
end


%% Ploting 
figure(2);
plot(Injection_range, frequency)
xlabel('I_{0} [mA/mm^{2}]')
ylabel('Firing Frquency [Hz]')
axis([.040 .1 0 .080])

%end

function setGlobalx(Iinj,duration,time)
global Iinj1 duration1 time1
Iinj1 = Iinj;
duration1 = duration;
time1 = time;

end
