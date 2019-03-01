% Author: Prabhu Chandhar, Chandhar Research Labs, Chennai, India.

clc;
clear all;

format long;
f = 2.5e9; % Hz
c = 3e8; % m/s
lambda = c/f; % m
del = .5*lambda;


%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 1a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = 10;
M =  100;
rhoudB = 10;
rhou = 10.^(rhoudB/10);
rhopdB = [-10:3:30];
lp = 1;
for  rhop = 10.^(rhopdB/10)
    pu_p = rhop*(4*pi*500/lambda)^2;
    [Exp_SINR_exact(lp),Exp_SINR_approx(lp),Exp_InvSINR_exact(lp),Exp_InvSINR_approx(lp),ErgodicRate_Simulation_ZF(lp),ErgodicRate_Simulation_MRC(lp)] = Calc_SPAWC2018(M,K,del,lambda,pu_p,rhou);
    lp = lp+1;
end

R_ub_exact = log2(1+Exp_SINR_exact);
R_ub_approx = log2(1+Exp_SINR_approx);
R_Simulation = ErgodicRate_Simulation_ZF;
R_lb_approx = log2(1+1./Exp_InvSINR_approx);
R_lb_exact = log2(1+1./Exp_InvSINR_exact);


figure(1);plot(rhopdB,R_ub_exact,'r');hold on;
figure(1);plot(rhopdB,R_ub_approx,'r--');hold on;
figure(1);plot(rhopdB,R_Simulation,'m');hold on;
figure(1);plot(rhopdB,R_lb_exact,'b');hold on;
figure(1);plot(rhopdB,R_lb_approx,'b--');hold on;
xlabel('Pilot SNR [\rho_p] (dB)');ylabel('Ergodic rate (bits/sec/Hz)');
legend(5,'Upper bound (Exact-SINR)','Upper bound (Approximate-SINR)','Simulation','Lower bound (Exact-SINR)','Upper bound (Approximate-SINR)');

%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 1b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear Exp_SINR_exact Exp_SINR_approx Exp_InvSINR_exact Exp_InvSINR_approx ErgodicRate_Simulation_ZF ErgodicRate_Simulation_MRC;

M =  500;
rhoudB = 10;
rhou = 10.^(rhoudB/10);
rhopdB = [-30:3:30];
lp = 1;
for  rhop = 10.^(rhopdB/10)
    pu_p = rhop*(4*pi*500/lambda)^2;
    [Exp_SINR_exact(lp),Exp_SINR_approx(lp),Exp_InvSINR_exact(lp),Exp_InvSINR_approx(lp),ErgodicRate_Simulation_ZF(lp),ErgodicRate_Simulation_MRC(lp)] = Calc_SPAWC2018(M,K,del,lambda,pu_p,rhou);
    lp = lp+1;
end


R_ub_exact = log2(1+Exp_SINR_exact);
R_ub_approx = log2(1+Exp_SINR_approx);
R_Simulation = ErgodicRate_Simulation_ZF;
R_lb_approx = log2(1+1./Exp_InvSINR_approx);
R_lb_exact = log2(1+1./Exp_InvSINR_exact);


figure(2);plot(rhopdB,R_ub_exact,'r');hold on;
figure(2);plot(rhopdB,R_ub_approx,'r--');hold on;
figure(2);plot(rhopdB,R_Simulation,'m');hold on;
figure(2);plot(rhopdB,R_lb_exact,'b');hold on;
figure(2);plot(rhopdB,R_lb_approx,'b--');hold on;

xlabel('Pilot SNR [\rho_p] (dB)');ylabel('Ergodic rate (bits/sec/Hz)');
legend(5,'Upper bound (Exact-SINR)','Upper bound (Approximate-SINR)','Simulation','Lower bound (Exact-SINR)','Upper bound (Approximate-SINR)');

%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear Exp_SINR_exact Exp_SINR_approx Exp_InvSINR_exact Exp_InvSINR_approx ErgodicRate_Simulation_ZF ErgodicRate_Simulation_MRC;

M_array =  50:50:500;
rhou =  10;
rhopdB = 30;
rhop = 10.^(rhopdB/10);
pu_p = rhop*(4*pi*500/lambda)^2;

lp = 1;
for  M = [50:50:500]
    [Exp_SINR_exact(lp),Exp_SINR_approx(lp),Exp_InvSINR_exact(lp),Exp_InvSINR_approx(lp),ErgodicRate_Simulation_ZF(lp),ErgodicRate_Simulation_MRC(lp)] = Calc_SPAWC2018(M,K,del,lambda,pu_p,rhou);
    lp = lp+1;
end


R_ub_exact = log2(1+Exp_SINR_exact);
R_ub_approx = log2(1+Exp_SINR_approx);
R_Simulation = ErgodicRate_Simulation_ZF;
R_lb_approx = log2(1+1./Exp_InvSINR_approx);
R_lb_exact = log2(1+1./Exp_InvSINR_exact);


figure(3);plot(M_array,R_ub_exact,'r');hold on;
figure(3);plot(M_array,R_ub_approx,'r--');hold on;
figure(3);plot(M_array,R_Simulation,'m');hold on;
figure(3);plot(M_array,R_lb_exact,'b');hold on;
figure(3);plot(M_array,R_lb_approx,'b--');hold on;

xlabel('Number of GS antennas [M]');ylabel('Ergodic rate (bits/sec/Hz)');
legend(5,'Upper bound (Exact-SINR)','Upper bound (Approximate-SINR)','Simulation','Lower bound (Exact-SINR)','Upper bound (Approximate-SINR)');

%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear Exp_SINR_exact Exp_SINR_approx Exp_InvSINR_exact Exp_InvSINR_approx ErgodicRate_Simulation_ZF ErgodicRate_Simulation_MRC;

K = 10;
M =  100;
rhopdB = 30;
rhop = 10.^(rhopdB/10);
rhoudB = [-30:3:30];
lp = 1;
for  rhou = 10.^(rhoudB/10)
    pu_p = rhop*(4*pi*500/lambda)^2;
    [Exp_SINR_exact(lp),Exp_SINR_approx(lp),Exp_InvSINR_exact(lp),Exp_InvSINR_approx(lp),ErgodicRate_Simulation_ZF(lp),ErgodicRate_Simulation_MRC(lp)] = Calc_SPAWC2018(M,K,del,lambda,pu_p,rhou);
    lp = lp+1;
end
R_Simulation_ZF = ErgodicRate_Simulation_ZF;
R_Simulation_MRC = ErgodicRate_Simulation_MRC;

figure(4);plot(rhoudB,R_Simulation_ZF);hold on;
figure(4);plot(rhoudB,R_Simulation_MRC,'r');hold on;

rhopdB = -10;
rhop = 10.^(rhopdB/10);
rhoudB = [-30:3:30];
lp = 1;
for  rhou = 10.^(rhoudB/10)
    pu_p = rhop*(4*pi*500/lambda)^2;
    [Exp_SINR_exact(lp),Exp_SINR_approx(lp),Exp_InvSINR_exact(lp),Exp_InvSINR_approx(lp),ErgodicRate_Simulation_ZF(lp),ErgodicRate_Simulation_MRC(lp)] = Calc_SPAWC2018(M,K,del,lambda,pu_p,rhou);
    lp = lp+1;
end
R_Simulation_ZF = ErgodicRate_Simulation_ZF;
R_Simulation_MRC = ErgodicRate_Simulation_MRC;

figure(4);plot(rhoudB,R_Simulation_ZF,'b--');hold on;
figure(4);plot(rhoudB,R_Simulation_MRC,'r--');hold on;
xlabel('Data SNR [\rho_u] (dB)');ylabel('Ergodic rate (bits/sec/Hz)');
legend('ZF (\rho_p = 30 dB)','MRC (\rho_p = 30 dB)','ZF (\rho_p = -10 dB)','MRC (\rho_p = -10 dB)');
