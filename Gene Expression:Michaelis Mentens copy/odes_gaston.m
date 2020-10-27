function dydt = odes_gaston(t, y, params)

%% Unpack params

ktx = params(1);
%kd = params(2);
%n = params(3);
delta_mRNA = params(4);
%TF = params(5);
ktl = params(6);
delta_Protein = params(7);
mu = params(8);
kcat = params(9);
KiCP = params(10);
Korn = params(11);
KCP = params(12);
KiPi = params(13);
%delta_c = params(14);
%% Populate variables

%%% Mutation
mRNA = y(1);


%%% Building blocks
protein = y(2);


%%% Product
citrulline = y(3);

orn = y(4);

cp = y(5);

Pi = y(6);
%% ODEs 

dydt(size(y,1),1) = 0; % Setup ODE structure

% mRNA
dydt(1) = ktx - (delta_mRNA + mu) * mRNA; %mRNA conc. ODE
% protein  
dydt(2) = ktl * mRNA - (delta_Protein + mu) * protein ;    %Protein conc. ODE
%citrulline
dydt(3) = (kcat * protein * cp * orn) / (KiCP * Korn * (1 + Pi/KiPi) + KCP * orn * (1 + Pi/KiPi) + Korn * cp + cp * orn); %- delta_c * citrulline;

dydt(4) = -(kcat * protein * cp * orn) / (KiCP * Korn * (1 + Pi/KiPi) + KCP * orn * (1 + Pi/KiPi) + Korn * cp + cp * orn) ;

dydt(5) = -(kcat * protein * cp * orn) / (KiCP * Korn * (1 + Pi/KiPi) + KCP * orn * (1 + Pi/KiPi) + Korn * cp + cp * orn) ;


dydt(6) = (kcat * protein * cp * orn) / (KiCP * Korn * (1 + Pi/KiPi) + KCP * orn * (1 + Pi/KiPi) + Korn * cp + cp * orn);

 



