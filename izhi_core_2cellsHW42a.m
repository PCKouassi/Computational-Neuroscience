function izhi_core_2cellsHW42a()

% solves the Izhikeviclearch using forward Euler

global tau_1_exc tau_2_exc Pnorm_exc Ps_exc A_exc B_exc 
global tau_1_inh tau_2_inh Ps_inh  Pnorm_inh A_inh B_inh
global c d Vrest Vthresh delta_t tmax Y_Init Delta Pnorm tau1 tau2

set(0,'DefaultAxesFontSize',10,'defaultaxeslinewidth',2,...
    'defaultlinelinewidth',2.,'defaultpatchlinewidth',1.5)

%% System Parameters
% parameter sets in Iz03
pars=[0.02      0.2     -65      6       14 ;...    % tonic spiking
    0.02      0.25    -65      6       0.5 ;...   % phasic spiking
    0.02      0.2     -50      2       15 ;...    % tonic bursting
    0.02      0.25    -55     0.05     0.6 ;...   % phasic bursting
    0.02      0.2     -55     4        10 ;...    % mixed mode
    0.01      0.2     -65     8        30 ;...    % spike frequency adaptation
    0.02      -0.1    -55     6        0  ;...    % Class 1
    0.2       0.26    -65     0        0  ;...    % Class 2
    0.02      0.2     -65     6        7  ;...    % spike latency
    0.05      0.26    -60     0        0  ;...    % subthreshold oscillations
    0.1       0.26    -60     -1       0  ;...    % resonator Iinj=1.2 resonance at deltat=50
    0.02      -0.1    -55     6        0  ;...    % integrator
    0.03      0.25    -60     4        0;...      % rebound spike
    0.03      0.25    -52     0        0;...      % rebound burst
    0.03      0.25    -60     4        0  ;...    % threshold variability
    1         1.5     -60     0      -65  ;...    % bistability
    1       0.2     -60     -21      0  ;...    % DAP
    0.02      1       -55     4        0  ;...    % accomodation
    -0.02      -1      -60     8        80 ;...    % inhibition-induced spiking
    -0.026     -1      -45     0        80];       % inhibition-induced bursting


%% Synaptic Parameters

trise=tau1*tau2/(tau1-tau2);
Pnorm=((tau2/tau1)^(trise/tau1)-(tau2/tau1)^(trise/tau2))^-1;

trise_exc=tau_1_exc*tau_2_exc/(tau_1_exc-tau_2_exc);
Pnorm_exc=((tau_2_exc/tau_1_exc)^(trise_exc/tau_1_exc)-(tau_2_exc/tau_1_exc)^(trise_exc/tau_2_exc))^-1;

trise_inh=tau_1_inh*tau_2_inh/(tau_1_inh-tau_2_inh);
Pnorm_inh=((tau_2_inh/tau_1_inh)^(trise_inh/tau_1_inh)-(tau_2_inh/tau_1_inh)^(trise_inh/tau_2_inh))^-1;

%% pre-allocate memory
maxnumsteps=ceil(tmax/delta_t+1);
t(1:maxnumsteps)=NaN;
Y(1:maxnumsteps,length(Y_Init))=NaN;
spikes1=NaN(1,1000);

% need integer ratio tmax/dt
nsteps=tmax/delta_t;
if (floor(nsteps)~=nsteps)
    fprintf(' not an integer number of steps nsteps = %g \n',nsteps);
end

Y(1,:)= Y_Init;

k=1;
t(1)=0;
%% Simulation
ispike1=0;

exc_num = [];
inh_num = [];
while t(k) < tmax
    [Y(k+1,:)]=FE(t(k),Y(k,:),delta_t);
    [exc_num(k), inh_num(k)]=increment(Ps_exc,Ps_inh);
    Y(k+1,5) = A_exc+Y(k+1, 5);
    Y(k+1,6) = B_exc+Y(k+1, 6);
    Y(k+1,7) = A_inh+Y(k+1, 7);
    Y(k+1,8) = B_inh+Y(k+1, 8);
    % resetting when neuron 1 spikes
    if Y(k+1,1)>=Vthresh
        ispike1=ispike1+1;
        spikes1(ispike1)=t(k);
        Y(k+1,1)=c; % Reset volatage of neuron to c value
        Ps1=Pnorm*(Y(k+1,3)-Y(k+1,4));
        Y(k+1,2)=Y(k+1,2)+d;
        Y(k+1,3)=Y(k+1,3)+Delta*(1-Ps1);
        Y(k+1,4)=Y(k+1,4)+Delta*(1-Ps1);
        
    end

    t(k+1)=t(k)+delta_t;
    k=k+1;

    if mod(k,10000) == 0
        fprintf('%g steps t = %g dt =%g \n',k,t(k),delta_t);
    end

end


%% Determine periods
V1=Y(:,1);
Vsignchange=V1(1:end-1).*V1(2:end);
indsignchange=find(Vsignchange<0);
periods=indsignchange(3:end)-indsignchange(1:end-2);
periods=delta_t*periods;

%% Data Analysis
tmini= 20;
for i = 1:tmax/tmini
   nspike(i) = length(find(spikes1> tmini*(i-1) & spikes1 < tmini*(i)));
   rateSection(i) = nspike(i)/tmini;
end

global mean_spike variance isi total_exc_r total_inh_r total_neuron_r
% Stats
mean_spike = mean(nspike);
variance = var(nspike);

% ISI
isi = spikes1(2:end)-spikes1(1:end-1);

% Firing Rate Calcs
total_nspike = sum(nspike);
total_m_exc = sum(exc_num);
total_m_inh = sum(inh_num);

total_neuron_r = total_nspike/tmax;
total_exc_r  =  total_m_exc/tmax;
total_inh_r  =  total_m_inh/tmax;

%% Plotting Parameters
set(0,'DefaultAxesFontSize',15,'defaultaxeslinewidth',2,...
    'defaultlinelinewidth',2.,'defaultpatchlinewidth',1.5)

%% Plotting

figure(1);
subplot(2,1,1)
plot(t(1:nsteps),Y(1:nsteps,1),'b');
title('Evolution of Membrane Potential')
xlabel('Time [ms]')
ylabel('Membrane Potential [mV]')

diff_exc = Y(1:nsteps,5)-Y(1:nsteps,6);
diff_inh =  Y(1:nsteps,7)-Y(1:nsteps,8);
subplot(2,1,2)
plot(t(1:nsteps),diff_exc,'b',t(1:nsteps),diff_inh,'r');
title('Evolution of A-B ')
xlabel('Time [ms]')
ylabel('A-B Difference')
legend('Exc','Inh')


figure(2)
subplot(2,1,1)
histogram(isi,100)
set(gca,'YScale','linear')
title('ISI Histogram')
xlabel('ISI')
ylabel('Linear')

subplot(2,1,2)
histogram(isi,100)
set(gca,'YScale','log')
title('IS Histogram')
xlabel('ISI')
ylabel('Log')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Yp=F(time,Y)

global G_exc Esyn_exc Pnorm_exc Ps_exc  tau_1_exc tau_2_exc
global G_inh Esyn_inh Pnorm_inh Ps_inh  tau_1_inh tau_2_inh
global tau1 tau2 a b 

V_neuron = Y(1);
u_neuron = Y(2);
A_neuron = Y(3);
B_neuron = Y(4);
A_exc = Y(5);
B_exc = Y(6);
A_inh = Y(7);
B_inh = Y(8);

dA1_dt=-A_neuron/tau1;
dB1_dt=-B_neuron/tau2;

%Ps_neuron=Pnorm*(A1-B1);
Ps_exc = Pnorm_exc*(A_exc-B_exc);
Ps_inh = Pnorm_inh*(A_inh-B_inh);

%Calculate voltage of postynaptic neuorn.
Yp(1)=0.04*V_neuron^2+5*V_neuron+140-G_exc*Ps_exc*(V_neuron-Esyn_exc)-G_inh*Ps_inh*((V_neuron-Esyn_inh));

Yp(2)=a*(b*V_neuron-u_neuron);
Yp(3)=dA1_dt;
Yp(4)=dB1_dt;
Yp(5) = -A_exc/tau_1_exc;
Yp(6) = -B_exc/tau_2_exc;
Yp(7) = -A_inh/tau_1_inh;
Yp(8) = -B_inh/tau_2_inh;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y]=FE(time,Y,dt)

global model

Y=Y+dt*F(time,Y);

end

%% Increment Function
function[m_exc, m_inh] = increment(Ps_exc,Ps_inh)
global N_exc r_exc delta_Ps_exc
global N_inh r_inh delta_Ps_inh
global delta_t
global A_exc B_exc A_inh B_inh

m_exc = binornd(N_exc, (r_exc*delta_t));
m_inh = binornd(N_inh, (r_inh*delta_t));

A_exc = m_exc*delta_Ps_exc*(1-Ps_exc);
B_exc = m_exc*delta_Ps_exc*(1-Ps_exc);

A_inh = m_inh*delta_Ps_inh*(1-Ps_inh);
B_inh = m_inh*delta_Ps_inh*(1-Ps_inh);

end