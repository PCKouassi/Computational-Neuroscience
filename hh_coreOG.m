function hh_coreOG

% solves the Hodgkin-Huxley equations using forward Euler

global Ena Ek Eca El Gna Gk Gl Gca
global Iinj0 Iinj1 Iinj2 Iinj3 Iinj4 time1 time2 time3 time4 duration1 duration2 duration3  duration4

set(0,'defaultaxesfontsize',16,'defaultaxeslinewidth',1.5,...
    'defaultlinelinewidth',2.,'defaultpatchlinewidth',1.5)

dt=.01; t_min=0;  t_max=100;

method='ode45';
%method='FE';
fprintf(' Method is %s\n',method)

% baseline injected current
Iinj0=0;

% generate additional brief pulse of injected current
time1=5;
duration1=50;
Iinj1=0.03; % in mA/mm^2

time2=70;
duration2=5;
Iinj2=0.03; % in mA/mm^2

% time3=520;
% duration3=5;
% Iinj3=0.2; % in mA/mm^2

% time4=540;
% duration4=5;
% Iinj4=0.2; % in mA/mm^2


% rough approximation of the rest state

V(1)=-60.75; n(1)=0.384; m(1)=0.0863; h(1)=0.444; M(1)=.353; H(1)=0.00628;

% initial condition
Y(1,:)=[V(1),n(1),m(1),h(1),M(1),H(1)];


switch method
    case 'FE'
        tic
        [time,Y]=FE(@F_hh,[t_min,t_max],Y(1,:),dt);
        toc
    case 'ode45'
        tic
        tol=1e-11;
        opts = odeset('RelTol',tol,'AbsTol',tol,'MaxStep',duration1/2);
        [time,Y]=ode15s(@F_hh,[t_min,t_max],Y(1,:),opts);
        toc
end

V=Y(:,1); n_gate=Y(:,2); m_gate=Y(:,3); h_gate=Y(:,4); M_gate=Y(:,5); H_gate = Y(:,6);

% determine periods
Vsignchange=V(1:end-1).*V(2:end);
ind_sign_change=find(Vsignchange<=0);
ind_voltage_negative=find(V<0);
ind_rise=intersect(ind_sign_change,ind_voltage_negative);
spike_times=time(ind_rise)
if (length(spike_times)>1)
    periods=spike_times(2:end)-spike_times(1:end-1)
else
    fprintf(' Less than 2 spikes \n');
end

total_current = Gl.*(V-El)+Gk.*n_gate.^4.*(V-Ek)+Gna*m_gate.^3.*h_gate.*(V-Ena)

current_Na = Gna*m_gate.^3.*h_gate.*(V-Ena);
current_K = Gk.*n_gate.^4.*(V-Ek);
current_Ca = Gca.*M_gate.^2.*H_gate.*(V-Eca);
current_leak = Gl.*(V-El);

% plotting

figure(1);

subplot(3,1,1)
plot(time,V,'-');
axis([0 t_max -80 50]);
ylabel('Voltage [mV]')
xlabel('Time [ms]')
legend(' voltage');

subplot(3,1,2)
plot(time,n_gate,'m',time,m_gate,'r',time,h_gate,'b');
xlabel('Time [ms]')
ylabel('Gate Value')
legend('n','m','h');
axis([0 t_max 0 1]);

subplot(3,1,3)
plot(time,current_Na,'m', time,current_K,'r',time,current_leak,'b',time, total_current, 'k');
legend('INa','IK','IL','Total Current');
ylabel('Current [mA]')
xlabel('Time [ms]')
axis([0 t_max 0 1]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time,Y]=FE(RHS,time_interval,Y_init,dt)

t_min=time_interval(1); t_max=time_interval(2);

% pre-allocate memory for speed
maxnumsteps=ceil((t_max-t_min)/dt)+1;
time=NaN(1,maxnumsteps); Y=NaN(maxnumsteps,4);

k=1; time(k)=t_min; Y(k,:)=Y_init;

while time(k) < t_max
    
    if (time(k)>t_max-dt+1e-10)
        dt=t_max-time(k);
        fprintf(' adjusting the last time step dt = %16.14g time = %16.14g \n',dt,time(k))
    end
    [Y(k+1,:)]=Y(k,:)+dt*RHS(time(k),Y(k,:))';
    time(k+1)=time(k)+dt;
    k=k+1;
    
    if mod(k,50000) == 0
        fprintf('%g steps t = %g V =%g \n',k,time(k),Y(k,1));
    end
    
end

% in case the preallocation allocated additional entries, remove them
time=time(1:k-1)'; Y=Y(1:k-1,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Yp = F_hh(time,Y)

global Ena Ek Eca El Gna Gk Gl Gca
global Iinj0 Iinj1 Iinj2 Iinj3 Iinj4 time1 time2 time3 time4 duration1 duration2 duration3  duration4

V=Y(1);
n=Y(2);
m=Y(3);
h=Y(4);
M=Y(5);
H=Y(6);

% injected current
Iinj=Iinj0+ Iinj1*(time>=time1)*(time<=time1+duration1)+Iinj2*(time>=time2)*(time<=time2+duration2)%+Iinj3*(time>=time3)*(time<=time3+duration3)%+Iinj4*(time>=time4)*(time<=time4+duration4);

Gna=1.20; % in mS/mm^2
%Gna=0;
Ena=50;

Gk=0.36;%1.5;%0.36;
%Gk=0;
Ek=-77;

Gca = 0.015;
Eca = 120; %mV;

Cm=0.01; % micro Farad/mm^2; time therefore measured in ms
Gl=0.003;
El=-54.387;

Vthresh=-40;
alpham=0.1*(V-Vthresh)/(1-exp(-0.1*(V-Vthresh)));
betam=4*exp(-0.0556*(V+65));
mp=alpham*(1-m)-betam*m;

alphah=0.07*exp(-0.05*(V+65));
betah=1/(1+exp(-0.1*(V+35)));
hinfinity=alphah/(alphah+betah);
tauh=1/(alphah+betah);
hp=(hinfinity-h)/tauh;

Vthresh=-55;
alphan=0.01*(V-Vthresh)/(1-exp(-0.1*(V-Vthresh)));
betan=0.125*exp(-0.0125*(V+65));
ninfinity=alphan/(alphan+betan);
taun=1/(alphan+betan);
np=(ninfinity-n)/taun;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Minfinity= 1/(1+exp(-((V+57)/6.2)));
tauM = 1/(exp(-((V+132)/16.7))+exp(((V+16.8)/18.2)));
Mp = Minfinity-M/tauM;

Hinfinity= 1/(1+exp((V+81)/4));
if V<-80
    tauH = exp((V+467)/66.6);
elseif V>= -80
    tauH = exp(-(V+22)/10.5)+28;
end
Hp= (Hinfinity-H)/tauH;

Vp=1/Cm*(Gna*m^3*h*(Ena-V)+Gk*n^4*(Ek-V)+Gl*(El-V))

Yp(1)=Vp;
Yp(2)=np;
Yp(3)=mp;
Yp(4)=hp;
Yp(5)= Mp;
Yp(6)= Hp;
Yp=Yp';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
