%{
Meeus Ephemeris funciton. Takes in the time in JDE and computes the orbital
elements for the planets at that time

Inputs: Time in JDE (days since J2000 epoch + 2451545)
Outputs: Orbital elements for the 8 planets plus Pluto
%}
function [Ephem] = Meeus(JDE)
AU = 1.49597870700e8;

T = (JDE - 2451545.0)/36525;
%% Mercury Ephemeris data
Ephem.Mercury.L = 252.250906 + (149472.6746358)*T + (-0.00000535)*T^2 + (0.000000002)*T^3;
Ephem.Mercury.a = 0.387098310*AU;
Ephem.Mercury.e = 0.20563175 + (0.000020406)*T + (-0.0000000284)*T^2 + (-0.00000000017)*T^3;
Ephem.Mercury.i = 7.004986 + (-0.0059516)*T + (0.00000081)*T^2 + (0.000000041)*T^3;
Ephem.Mercury.Omega = wrapTo360(48.330893 + (-0.1254229)*T + (-0.00008833)*T^2 + (-0.000000196)*T^3);
Ephem.Mercury.PI = wrapTo360(77.456119 + (0.1588643)*T + (-0.00001343)*T^2 + (0.000000039)*T^3);
Ephem.Mercury.omega = wrapTo360(Ephem.Mercury.PI - Ephem.Mercury.Omega);
Ephem.Mercury.M = wrapTo360(Ephem.Mercury.L - Ephem.Mercury.PI);
C_cen = (2*(Ephem.Mercury.e) - ((Ephem.Mercury.e)^3)/4 + (5/96)*(Ephem.Mercury.e)^5)*sind(Ephem.Mercury.M) + ((5/4)*(Ephem.Mercury.e)^2 - (11/24)*(Ephem.Mercury.e)^4)*sind(2*Ephem.Mercury.M) + ...
    ((13/12)*(Ephem.Mercury.e)^3 - (43/64)*(Ephem.Mercury.e)^5)*sind(3*Ephem.Mercury.M) + (103/96)*(Ephem.Mercury.e)^4*sind(4*Ephem.Mercury.M) + (1097/960)*(Ephem.Mercury.e)^5*sind(5*Ephem.Mercury.M);
Ephem.Mercury.nu = wrapTo360(Ephem.Mercury.M + (180/pi)*C_cen);
Ephem.Mercury.w = wrapTo360(Ephem.Mercury.PI - Ephem.Mercury.Omega);
%% Venus Ephemeris data
Ephem.Venus.L = 181.979801 + (58517.8156760)*T + (0.00000165)*T^2 + (-0.000000002)*T^3;
Ephem.Venus.a = 0.72332982*AU;
Ephem.Venus.e = 0.00677188 + (-0.000047766)*T + (0.0000000975)*T^2 + (0.00000000044)*T^3;
Ephem.Venus.i = 3.394662 + (-0.0008568)*T + (-0.00003244)*T^2 + (0.000000010)*T^3;
Ephem.Venus.Omega = wrapTo360(76.679920 + (-0.2780080)*T + (-0.00014256)*T^2 + (-0.000000198)*T^3);
Ephem.Venus.PI = wrapTo360(131.563707 + (0.0048646)*T + (-0.00138232)*T^2 + (-0.000005332)*T^3);
Ephem.Venus.M = wrapTo360(Ephem.Venus.L - Ephem.Venus.PI);
C_cen = (2*(Ephem.Venus.e) - ((Ephem.Venus.e)^3)/4 + (5/96)*(Ephem.Venus.e)^5)*sind(Ephem.Venus.M) + ((5/4)*(Ephem.Venus.e)^2 - (11/24)*(Ephem.Venus.e)^4)*sind(2*Ephem.Venus.M) + ...
    ((13/12)*(Ephem.Venus.e)^3 - (43/64)*(Ephem.Venus.e)^5)*sind(3*Ephem.Venus.M) + (103/96)*(Ephem.Venus.e)^4*sind(4*Ephem.Venus.M) + (1097/960)*(Ephem.Venus.e)^5*sind(5*Ephem.Venus.M);
Ephem.Venus.nu = wrapTo360(Ephem.Venus.M + (180/pi)*C_cen);
Ephem.Venus.w = wrapTo360(Ephem.Venus.PI - Ephem.Venus.Omega);
%% Earth Ephemeris data
Ephem.Earth.L = 100.466449 + (35999.3728519)*T + (-0.00000568)*T^2 + (0.0)*T^3;
Ephem.Earth.a = 1.000001018*AU;
Ephem.Earth.e = 0.01670862 + (-0.000042037)*T + (-0.0000001236)*T^2 + (0.00000000004)*T^3;
Ephem.Earth.i = 0 + (0.0130546)*T + (-0.00000931)*T^2 + (-0.000000034)*T^3;
Ephem.Earth.Omega = wrapTo360(174.873174 + (-0.2410908)*T + (0.00004067)*T^2 + (-0.000001327)*T^3);
Ephem.Earth.PI = wrapTo360(102.937348 + (0.3225557)*T + (0.00015026)*T^2 + (0.000000478)*T^3);
Ephem.Earth.M = wrapTo360(Ephem.Earth.L - Ephem.Earth.PI);
C_cen = (2*(Ephem.Earth.e) - ((Ephem.Earth.e)^3)/4 + (5/96)*(Ephem.Earth.e)^5)*sind(Ephem.Earth.M) + ((5/4)*(Ephem.Earth.e)^2 - (11/24)*(Ephem.Earth.e)^4)*sind(2*Ephem.Earth.M) + ...
    ((13/12)*(Ephem.Earth.e)^3 - (43/64)*(Ephem.Earth.e)^5)*sind(3*Ephem.Earth.M) + (103/96)*(Ephem.Earth.e)^4*sind(4*Ephem.Earth.M) + (1097/960)*(Ephem.Earth.e)^5*sind(5*Ephem.Earth.M);
Ephem.Earth.nu = wrapTo360(Ephem.Earth.M + (180/pi)*C_cen);
Ephem.Earth.w = wrapTo360(Ephem.Earth.PI - Ephem.Earth.Omega);
%% Mars Ephemeris data
Ephem.Mars.L = 355.433275 + (19140.2993313)*T + (0.00000261)*T^2 + (-0.000000003)*T^3;
Ephem.Mars.a = 1.523679342*AU;
Ephem.Mars.e = 0.09340062 + (0.000090483)*T + (-0.0000000806)*T^2 + (-0.00000000035)*T^3;
Ephem.Mars.i = 1.849726 + (-0.0081479)*T + (-0.00002255)*T^2 + (-0.000000027)*T^3;
Ephem.Mars.Omega = wrapTo360(49.558093 + (-0.2949846)*T + (-0.00063993)*T^2 + (-0.000002143)*T^3);
Ephem.Mars.PI = wrapTo360(336.060234 + (0.4438898)*T + (-0.00017321)*T^2 + (0.000000300)*T^3);
Ephem.Mars.M = wrapTo360(Ephem.Mars.L - Ephem.Mars.PI);
C_cen = (2*(Ephem.Mars.e) - ((Ephem.Mars.e)^3)/4 + (5/96)*(Ephem.Mars.e)^5)*sind(Ephem.Mars.M) + ((5/4)*(Ephem.Mars.e)^2 - (11/24)*(Ephem.Mars.e)^4)*sind(2*Ephem.Mars.M) + ...
    ((13/12)*(Ephem.Mars.e)^3 - (43/64)*(Ephem.Mars.e)^5)*sind(3*Ephem.Mars.M) + (103/96)*(Ephem.Mars.e)^4*sind(4*Ephem.Mars.M) + (1097/960)*(Ephem.Mars.e)^5*sind(5*Ephem.Mars.M);
Ephem.Mars.nu = wrapTo360(Ephem.Mars.M + (180/pi)*C_cen);
Ephem.Mars.w = wrapTo360(Ephem.Mars.PI - Ephem.Mars.Omega);
%% Jupiter Ephemeris data
Ephem.Jupiter.L = 34.351484 + (3034.9056746)*T + (-0.00008501)*T^2 + (0.000000004)*T^3;
Ephem.Jupiter.a = (5.202603191 + (0.0000001913)*T)*AU;
Ephem.Jupiter.e = 0.04849485 + (0.000163244)*T + (-0.0000004719)*T^2 + (-0.00000000197)*T^3;
Ephem.Jupiter.i = 1.303270 + (-0.0019872)*T + (0.00003318)*T^2 + (0.000000092)*T^3;
Ephem.Jupiter.Omega = wrapTo360(100.464441 + (0.1766828)*T + (0.00090387)*T^2 + (-0.000007032)*T^3);
Ephem.Jupiter.PI = wrapTo360(14.331309 + (0.2155525)*T + (0.00072252)*T^2 + (-0.000004590)*T^3);
Ephem.Jupiter.M = wrapTo360(Ephem.Jupiter.L - Ephem.Jupiter.PI);
C_cen = (2*(Ephem.Jupiter.e) - ((Ephem.Jupiter.e)^3)/4 + (5/96)*(Ephem.Jupiter.e)^5)*sind(Ephem.Jupiter.M) + ((5/4)*(Ephem.Jupiter.e)^2 - (11/24)*(Ephem.Jupiter.e)^4)*sind(2*Ephem.Jupiter.M) + ...
    ((13/12)*(Ephem.Jupiter.e)^3 - (43/64)*(Ephem.Jupiter.e)^5)*sind(3*Ephem.Jupiter.M) + (103/96)*(Ephem.Jupiter.e)^4*sind(4*Ephem.Jupiter.M) + (1097/960)*(Ephem.Jupiter.e)^5*sind(5*Ephem.Jupiter.M);
Ephem.Jupiter.nu = wrapTo360(Ephem.Jupiter.M + (180/pi)*C_cen);
Ephem.Jupiter.w = wrapTo360(Ephem.Jupiter.PI - Ephem.Jupiter.Omega);
%% Saturn Ephemeris data
Ephem.Saturn.L = 50.077471 + (1222.1137943)*T + (0.00021004)*T^2 + (-0.000000019)*T^3;
Ephem.Saturn.a = (9.554909596 + (-0.0000021389)*T)*AU;
Ephem.Saturn.e = 0.05550862 + (-0.000346818)*T + (-0.0000006456)*T^2 + (0.00000000338)*T^3;
Ephem.Saturn.i = 2.488878 + (0.0025515)*T + (-0.00004903)*T^2 + (0.000000018)*T^3;
Ephem.Saturn.Omega = wrapTo360(113.665524 + (-0.2566649)*T + (-0.00018345)*T^2 + (0.000000357)*T^3);
Ephem.Saturn.PI = wrapTo360(93.056787 + (0.5665496)*T + (0.00052809)*T^2 + (0.000004882)*T^3);
Ephem.Saturn.M = wrapTo360(Ephem.Saturn.L - Ephem.Saturn.PI);
C_cen = (2*(Ephem.Saturn.e) - ((Ephem.Saturn.e)^3)/4 + (5/96)*(Ephem.Saturn.e)^5)*sind(Ephem.Saturn.M) + ((5/4)*(Ephem.Saturn.e)^2 - (11/24)*(Ephem.Saturn.e)^4)*sind(2*Ephem.Saturn.M) + ...
    ((13/12)*(Ephem.Saturn.e)^3 - (43/64)*(Ephem.Saturn.e)^5)*sind(3*Ephem.Saturn.M) + (103/96)*(Ephem.Saturn.e)^4*sind(4*Ephem.Saturn.M) + (1097/960)*(Ephem.Saturn.e)^5*sind(5*Ephem.Saturn.M);
Ephem.Saturn.nu = wrapTo360(Ephem.Saturn.M + (180/pi)*C_cen);
Ephem.Saturn.w = wrapTo360(Ephem.Saturn.PI - Ephem.Saturn.Omega);
%% Pluto Ephemeris data
Ephem.Pluto.L = 238.92903833 + (145.20780515)*T;
Ephem.Pluto.a = AU*(39.48211675 + (-0.00031596)*T);
Ephem.Pluto.e = 0.24882730 + (0.00005170)*T;
Ephem.Pluto.i = 17.14001206 + (0.00004818)*T;
Ephem.Pluto.Omega = wrapTo360(110.30393684 + (-0.01183482)*T);
Ephem.Pluto.PI = wrapTo360(224.06891629 + (-0.04062942)*T);
Ephem.Pluto.M = wrapTo360(Ephem.Pluto.L - Ephem.Pluto.PI);
C_cen = (2*(Ephem.Pluto.e) - ((Ephem.Pluto.e)^3)/4 + (5/96)*(Ephem.Pluto.e)^5)*sind(Ephem.Pluto.M) + ((5/4)*(Ephem.Pluto.e)^2 - (11/24)*(Ephem.Pluto.e)^4)*sind(2*Ephem.Pluto.M) + ...
    ((13/12)*(Ephem.Pluto.e)^3 - (43/64)*(Ephem.Pluto.e)^5)*sind(3*Ephem.Pluto.M) + (103/96)*(Ephem.Pluto.e)^4*sind(4*Ephem.Pluto.M) + (1097/960)*(Ephem.Pluto.e)^5*sind(5*Ephem.Pluto.M);
Ephem.Pluto.nu = wrapTo360(Ephem.Pluto.M + (180/pi)*C_cen);
Ephem.Pluto.w = wrapTo360(Ephem.Pluto.PI - Ephem.Pluto.Omega);
end

