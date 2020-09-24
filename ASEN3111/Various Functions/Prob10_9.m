%%%%%%
%
%Conner Sahver
%ASEN 3111
%10.9
%Using .m files provided for equations
%
%%%%%%
close all
clear all
clc

%% c
%% Declaring Variables and givens
gamma = 1.4;
A2oAt = linspace(1, 1.4, 1000);
AeoAt = 1.53;
p01 = 1;
PE = 0.75;

%% Interate through shock locations and calculate the exit pressure and mach number for each one
for i = 1:length(A2oAt)
    [Msub(i) Msup(i)] = nozzle(A2oAt(i));
    [M2(i),p2op1(i),rho2orho1,t2ot1,deltasoR,p02op01(i)] = shock_calc(Msup(i));
    A2oAstar(i) = AratFindA(gamma, M2(i));
    AeoA2star(i) = AeoAt*(1/A2oAt(i))*A2oAstar(i); 
    [Msube(i) Msupe(i)] = nozzle(AeoA2star(i));
    [p02ope(i), t0ot, rho0orho] = isentropic(Msube(i));
    pe(i) = (1/p02ope(i))*p02op01(i)*p01;
end

%% Find the exit mach  number for the correct shock location
a = find(pe<=PE, 1, 'first');
Me = Msube(a)

%% d
clear all;
%% Declare variables and givens
gamma = 1.4;
A2oAt = linspace(1, 1.53, 1000);
AeoAt = 1.53;
p01 = 1;
PE = 0.154;

%% Interate through shock locations and calculate the exit pressure and mach number for each one
for i = 1:length(A2oAt)
    [Msub(i) Msup(i)] = nozzle(A2oAt(i));
    [M2(i),p2op1(i),rho2orho1,t2ot1,deltasoR,p02op01(i)] = shock_calc(Msup(i));
    A2oAstar(i) = AratFindA(gamma, M2(i));
    AeoA2star(i) = AeoAt*(1/A2oAt(i))*A2oAstar(i); 
    [Msube(i) Msupe(i)] = nozzle(AeoA2star(i));
    [p02ope(i), t0ot, rho0orho] = isentropic(Msupe(i));
    pe(i) = (1/p02ope(i))*p02op01(i)*p01;
end

%% Find the exit mach  number for the correct shock location
a = find(pe>=PE, 1, 'first');
Me = Msupe(a)