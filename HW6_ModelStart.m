function HW6_ModelStart

clear
clc
clf

set(0,'DefaultFigureWindowStyle','docked')
%% Initialize global variables
global modelSwitch niter dt tauw wmax wmin tautheta um rand

set(0,'DefaultAxesFontSize',15,'defaultaxeslinewidth',2,...
'defaultlinelinewidth',2.,'defaultpatchlinewidth',1.5)


%% Stimuli
u1 = [1.2,1,1]';
u2 = [1,1.5,1]';
u3 = [1,1.5,2]';
um = [u1 u2 u3];

%% Covariance & Eigan Value Calculations
covariance = cov(um');
[eigvector, eiganvalue]= eig(covariance);

%% Choose Model
modelSwitch = 0; % Value of 0 solves Covariance model, Value of 1 solves BCM model

%% Set parameters
tauw = 100000;
tautheta = 10000;
wmax = 200;
wmin = -200;
niter = 2E6
dt = 1;
%% Set initial values for variables w, v, Initialize Y
% wi = [0.01 0.01 0.01]
wi = eigvector(:,3)'
% wi = 10*eigvector(:,3)'
% wi = 200*eigvector(:,3)'

patterni = 3;
k = 1;

if modelSwitch == 0
    Y = zeros(niter, 5);
    vthresh = 0;
    Y(1,:) = [patterni, wi, vthresh];
else
    Y = zeros(niter, 5);
    vtheta_thresh = 0;
    Y(1,:) = [patterni, wi, vtheta_thresh];
end

%% Main
while k < niter
    [New_Y]= FE(Y(k,:),dt,k);

    New_Y(1) = rand;
    [Y(k+1,:)]= New_Y;
    if modelSwitch == 0
        if Y(k+1,2) > wmax
            Y(k+1,2) = wmax;
        elseif Y(k+1,2) < wmin
            Y(k+1,2) = wmin;
        end

        if Y(k+1,3) > wmax
            Y(k+1,3) = wmax;
        elseif Y(k+1,3) < wmin
            Y(k+1,3) = wmin;
        end

        if Y(k+1,4) > wmax
            Y(k+1,4) = wmax;
        elseif Y(k+1,4) < wmin
            Y(k+1,4) = wmin;
        end
    end
    k = k+1;
end

for i = 1:niter
    v_plot(i) = dot(Y(i,[2 3 4]), um(:,Y(i,1)));
end



%% Normalizng Weight Vectors
%1
normw1 = Y(:,2)/norm(Y(:,2));
normw2 = Y(:,3)/norm(Y(:,3));
normw3 = Y(:,4)/norm(Y(:,4));

normwfinal = round(Y(niter,[2 3 4])/norm(Y(niter,[2 3 4])), 2);

wi = round(wi,2);
%% Check if weights are different to inital condition

v1 = dot(normwfinal,u1);
v2 = dot(normwfinal,u2);
v3 = dot(normwfinal,u3);

final_rates = [v1, v2, v3];

fprintf('Final firirng rate = %f, %f, %f\n', final_rates(1), final_rates(2),final_rates(3))
fprintf('initial = %g \n',wi);
fprintf('final = %g \n',normwfinal);


%% Plotting
time = 1:niter;

figure(1)
subplot(4,1,1)
scatter(time, v_plot)
xlabel('Steps')
title('Firing Rate')

subplot(4,1,2)
plot(time,Y(:,2),time,Y(:,3),time,Y(:,4))
title('Weights')
xlabel('Steps')

subplot(4,1,3)
plot(time,normw1, time, normw2, time, normw3)
title('Normalized Weights')
xlabel('Steps')

subplot(4,1,4)
scatter(time, Y(:,5))
xlabel('Steps')
title('Threshold')

end

function Yp =F(Y,k)
global modelSwitch tauw tautheta um u rand
%% Differential Equations
% This section contains the differential equation(s) used to calculate the
% change in synaptic weights (w vector) at each time step

switch modelSwitch 
    case 0
        % Covariance Model
        rand = randi(3);
        u = um(:,rand);
        w = [Y(2),Y(3),Y(4)];
        v = dot(w',u);

        vthresh = (1/3)*(dot(w,um(:,1)') + dot(w,um(:,2)') + dot(w,um(:,3)'));
        new_weights = ((v-vthresh)*u)/tauw;
        Yp(1) = 0 ;
        Yp(2) = new_weights(1);
        Yp(3) = new_weights(2);
        Yp(4) = new_weights(3);
        Yp(5) = vthresh;

    case 1
        % BCM Model
        rand = randi(3);
        u = um(:,rand);

        w = [Y(2),Y(3),Y(4)];
        v = w*u;

        new_thresh = ((v^2)-Y(5))/tautheta;


        new_weights = (v*(v-Y(5)).*u)./tauw;

        Yp(1) = 0;
        Yp(2) = new_weights(1);
        Yp(3) = new_weights(2);
        Yp(4) = new_weights(3);
        Yp(5) = new_thresh;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end


function [Y]= FE(Y,dt,k)
Y=Y+dt*F(Y,k);
end




