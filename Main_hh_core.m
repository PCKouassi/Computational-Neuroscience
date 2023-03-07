function [periods] = Main_hh_core()

%%%%%%%%%%%%%%%%%% Global Function %%%%%%%%%%%%%%%%%%%%%%
% solves the Hodgkin-Huxley equations using forward Euler

global Ena Ek El Gna Gk Gl
global Iinj0 Iinj1 duration1 time1


set(0,'defaultaxesfontsize',16,'defaultaxeslinewidth',1.5,...
    'defaultlinelinewidth',2.,'defaultpatchlinewidth',1.5)


dt=.01; t_min=0;  t_max=800;

% Only using ode15s bc of stability
method='ode15s';
%method='FE';
fprintf(' Method is %s\n',method)

% baseline injected current
Iinj0=0;

% rough approximation of the rest state
V(1)=-65; n(1)=0.3; m(1)=0.1; h(1)=0.60;

% initial conditionc
Y(1,:)=[V(1),n(1),m(1),h(1)];

switch method
    case 'FE'
        tic
        [time,Y]=FE(@F_hh,[t_min,t_max],Y(1,:),dt);
        toc
    case 'ode15s'
        tic
        tol=1e-11;
        opts = odeset('RelTol',tol,'AbsTol',tol,'MaxStep',duration1/2);
        [time,Y]=ode15s(@F_hh,[t_min,t_max],Y(1,:),opts);
        toc
end

V=Y(:,1); n_gate=Y(:,2); m_gate=Y(:,3); h_gate=Y(:,4);

% determine periods
Vsignchange=V(1:end-1).*V(2:end);
ind_sign_change=find(Vsignchange<=0);
ind_voltage_negative=find(V<0);
ind_rise=intersect(ind_sign_change,ind_voltage_negative);
spike_times = time(ind_rise);

if (length(spike_times)>1)
    periods = spike_times(2:end)-spike_times(1:end-1);
else
    periods = [];
    fprintf(' Less than 2 spikes \n');
end

total_spike_time= sum(periods);

total_current = Gl.*(V-El)+Gk.*n_gate.^4.*(V-Ek)+Gna*m_gate.^3.*h_gate.*(V-Ena);

current_Na = Gna*m_gate.^3.*h_gate.*(V-Ena);
current_K = Gk.*n_gate.^4.*(V-Ek);
current_leak = Gl.*(V-El);

% plotting

% figure(1);
% 
% subplot(3,1,1)
% plot(time,V,'-');
% axis([0 t_max -80 50]);
% legend(' voltage');
% 
% subplot(3,1,2)
% plot(time,n_gate,'m',time,m_gate,'r',time,h_gate,'b');
% legend('n','m','h');
% axis([0 t_max 0 1]);
% 
% subplot(3,1,3)
% plot(time,current_Na,'m', time,current_K,'r',time,current_leak,'b',time, total_current, 'k');
% legend('INa','IK','IL','Total Current');
% axis([0 t_max 0 1]);

end



%%%%%%%%%%% FUNCTION DEFINITIONS %%%%%%%%%%%%%%%
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

global Ena Ek El Gna Gk Gl
global Iinj0 Iinj1 duration1 time1

V=Y(1);
n=Y(2);
m=Y(3);
h=Y(4);

% injected current
Iinj= Iinj0+ Iinj1*(time>=time1)*(time<=time1+duration1);

Gna=1.20; % in mS/mm^2
%Gna=0;
Ena=50;

Gk=0.36;%1.5;%0.36;
%Gk=0;
Ek=-77;

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

Vp=1/Cm*(Gna*m^3*h*(Ena-V)+Gk*n^4*(Ek-V)+Gl*(El-V)+Iinj);

Yp(1)=Vp;
Yp(2)=np;
Yp(3)=mp;
Yp(4)=hp;
Yp=Yp';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
