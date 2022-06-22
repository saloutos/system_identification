% script for processing arm log csv file for running the system ID optimizations
% Andrew SaLoutos, 6/22/2022

%% reset workspace
clear all
clc

% import data file
load("armlog2.mat");
data = armlog2;
% set output file name
output = "data/example_arm_logs/armlog2_processed_data.mat";

%% pull out q,qd,tau
t = data.t;
q = [data.q_BASE, data.q_SP, data.q_SR, data.q_ELBOW, data.q_WP, data.q_WR];
qd  = [data.qd_BASE, data.qd_SP, data.qd_SR, data.qd_ELBOW, data.qd_WP, data.qd_WR];
tau = [data.tau_BASE, data.tau_SP, data.tau_SR, data.tau_ELBOW, data.tau_WR, data.tau_WR1]; % two WR columns, first is actually WP and second is actually WR


%% plot everything
figure; plot(q(:,4)); % why is plotting being weird?

%% snip data
start_ind = 2700;
stop_ind = 10000;

t = t(start_ind:stop_ind)-t(start_ind);
q = q(start_ind:stop_ind,:);
qd = qd(start_ind:stop_ind,:);
tau = tau(start_ind:stop_ind,:);


%% filter q,qd,tau

% y = sgolayfilt(x,order,framelen,weights,dim)
filt_order = 5; % for cheetah leg, used order 4 and window length 101
window_size = 35;

qfilt = sgolayfilt(q, filt_order, window_size);
qdfilt = sgolayfilt(qd, filt_order, window_size);
taufilt = sgolayfilt(tau, filt_order, window_size);

%% finite difference to calculate qdd, filter qdd
dt = mean(diff(t));
qdd = gradient(qdfilt,dt);

%% save data
q = qfilt;
qd = qdfilt;
save(output, 't','q','qd','qdd','tau','taufilt');







