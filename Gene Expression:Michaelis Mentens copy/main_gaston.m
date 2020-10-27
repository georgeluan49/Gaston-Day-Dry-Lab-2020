%%Example code for Gaston High Scool
%% Written by Emma Albertini, Imperial College London
%% Edited and implemented by George Luan, Gaston Day School



%% Set-up to define style and color if you want your graphs to look nicer

close all; clear; clc;

set(0,'DefaultLineLineWidth',2);
set(0,'defaultAxesFontSize',24);
E_yellow = [0.88 0.78 0.02];
U_orange = [0.88 0.53 0];
M_red = [237, 12, 0]/255;
lightgrey = [0.85 0.85 0.85];
mediumgrey = [0.55 0.55 0.55];
darkgrey = [0.1 0.1 0.1];
R_blue = [0 114 189]/255;
G_green = [119 172 48]/255;
H_yellow = [237 177 32]/255;

tic



%% Parameters

%Transcription:
ktx = 2.203 * 10^(-7)  ;     %M/s transcription constant
kd = 0 ; %M  dissociation constant
n = 0  ;        %(no units) hill coefficient
delta_mRNA = 0.247; %/s degredation consant of mRNA
TF = 10^(-7);          %M concentration of transcirption factor

%Translation:
ktl = 0.2059 ;              %/s translation constant
delta_Protein = 0.008492;  %/s degredation constant for Protein

%Deluion
mu = 0.0072;

%citrulline:
kcat = 2940;
KiCP = 0.151 .* 10^(-3);
Korn = 1.1 .* 10^(-3);
KCP = 0.09 .* 10^(-3);
KiPi = 3.8 .* 10^(-3);
delta_c = 0.0116;
params = [ktx, kd, n, delta_mRNA, TF, ktl, delta_Protein, mu, kcat, KiCP, Korn, KCP, KiPi, delta_c];      





%% Initial variables
mRNA_0= 0;
protein_0=0;
citrulline_0 = 0;
cp_0 = (1.5) .* 10^(-4);
orn_0 = (4.0) .* 10^(-4);
Pi_0 = 0;


var_0 = [mRNA_0, protein_0, citrulline_0, cp_0, orn_0, Pi_0];

%% Populate variables
tspan = [0 3600];
options = odeset('NonNegative',1);
[t,y] = ode15s(@(t,y)odes_gaston(t, y, params),...
               tspan, ...
               var_0, options); 

%%% Mutation
mRNA = y(:,1);
protein = y(:,2);
citrulline = y(:,3);
cp = y(:,4);
orn = y(:,5);
Pi = y(:,6);
%%%% EM plot

close all

set(gcf, 'Position',  [100, 100, 800, 640])

plot(t, mRNA)

 hold on

 plot(t, protein)
 
 %hold on
 
 %plot(t, orn)
 
 %hold on
 
 %plot(t, cp)
 
 %hold on
 
 %plot(t, citrulline)

 %axis([0, 3600, 0, 5e-4]) 