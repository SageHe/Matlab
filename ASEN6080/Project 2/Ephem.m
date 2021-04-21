function [R V mu_p] = Ephem(JD,Planet,Frame)
%--------------------------------------------------
%
%  Center of reference frames are the Sun/Solar system barycenter
%
%  Inputs:
%       JD = Julian Day (double)
%       Planet = number from 1-9
%                   1 - Mercury
%                   2 - Venus
%                   3 - Earth
%                   4 - Mars
%                   5 - Jupiter
%                   6 - Saturn
%                   7 - Uranus
%                   8 - Neptune
%                   9 - Pluto
%       Frame = Desired Reference Frame (str)
%                   'EME2000'
%                   'EMO2000'
%
%  Outputs:
%       R = Planet's Position (km)
%       V = Planet's Velocity (km/s)
%       mu_p = Planet's gravitational constant (km^3/s^2)
%
%--------------------------------------------------

% These are in EMO2000

T = (JD - 2451545.0)/36525;
deg2rad = pi/180;
AU = 149597870.700;
mu_s = 132712440017.987;

[coeff] = Ephemeride_Coeff(Planet);

Tvec = [1; T; T*T; T*T*T];

% Mean longitude of Planet
L = coeff.L*Tvec*deg2rad;

% Semimajor axis of the orbit
a = coeff.a*Tvec*AU;

% Eccentricity of the orbit
e = coeff.e*Tvec;

% Inclination of the orbit
inc = coeff.i*Tvec*deg2rad;

% Longitude of the Ascending Node
W = coeff.W*Tvec*deg2rad;

% Longitude of the Perihelion
P = coeff.P*Tvec*deg2rad;

% Argument of perihelion
w = P - W;

% Mean anomaly of orbit
M = L - P;

% True anomaly of orbit
Ccen = (2*e - e^3/4 + 5/96*e^5)*sin(M) + (5/4*e^2 - 11/24*e^4)*sin(2*M) + ...
    (13/12*e^3 - 43/64*e^5)*sin(3*M) + 103/96*e^4*sin(4*M) + ...
    1097/960*e^5*sin(5*M);

nu = M + Ccen;

[R V] = COEstoRV(a,e,inc,W,w,nu,mu_s);

% Convert to EME2000 if necessary
if strcmpi(Frame,'EME2000')
    theta = 23.4393*pi/180;
    C = [1 0 0;
        0 cos(theta) -sin(theta);
        0 sin(theta) cos(theta)];
    R = C*R;
    V = C*V;
end

mu_p = coeff.mu_p;

end

function [coeff] = Ephemeride_Coeff(Planet)

if  Planet == 1;
    % Venus
    L = [252.250906     149472.6746358 -0.00000535      0.000000002];
    a = [0.387098310    0.0             0.0             0.0];
    e = [0.20563175     0.000020406    -0.0000000284   -0.00000000017];
    i = [7.004986      -0.0059516       0.00000081      0.000000041];
    W = [48.330893     -0.1254229      -0.00008833     -0.000000196];
    P = [77.456119	 	0.1588643      -0.00001343      0.000000039];
    mu_p = 2.20320804864179e4;
elseif  Planet == 2;
    % Venus
    L = [181.979801     58517.8156760   0.00000165     -0.000000002];
    a = [0.72332982     0.0             0.0             0.0];
    e = [0.00677188    -0.000047766     0.0000000975	0.00000000044];
    i = [3.394662      -0.0008568      -0.00003244      0.000000010];
    W = [76.679920     -0.2780080      -0.00014256     -0.000000198];
    P = [131.563707     0.0048646      -0.00138232     -0.000005332];
    mu_p = 3.2485859882646e5;
elseif Planet == 3;
    % Earth
    L = [100.466449     35999.3728519  -0.00000568      0.0];
    a = [1.000001018    0.0             0.0             0.0];
    e = [0.01670862    -0.000042037    -0.0000001236	0.00000000004];
    i = [0.0            0.0130546      -0.00000931     -0.000000034];
    W = [174.873174    -0.2410908       0.00004067     -0.000001327];
    P = [102.937348     0.3225557       0.00015026      0.000000478];
    mu_p = 3.98600432896939e5;
elseif Planet == 4;
    % Mars
    L = [355.433275     19140.2993313	0.00000261     -0.000000003];
    a = [1.523679342    0.0             0.0             0.0];
    e = [0.09340062     0.000090483    -0.0000000806   -0.00000000035];
    i = [1.849726      -0.0081479      -0.00002255     -0.000000027];
    W = [49.558093     -0.2949846      -0.00063993     -0.000002143];
    P = [336.060234     0.4438898      -0.00017321      0.000000300];
    mu_p = 4.28283142580671e4;
elseif Planet == 5;
    % Jupiter
    L = [34.351484      3034.9056746   -0.00008501      0.000000004];
    a = [5.202603191	0.0000001913    0.0             0.0];
    e = [0.04849485     0.000163244	   -0.0000004719   -0.00000000197];
    i = [1.303270      -0.0019872       0.00003318      0.000000092];
    W = [100.464441     0.1766828       0.00090387     -0.000007032];
    P = [14.331309      0.2155525       0.00072252     -0.000004590];
    mu_p = 1.26712767857796e8;
elseif Planet == 6;
    % Saturn
    L = [50.077471      1222.1137943    0.00021004     -0.000000019];
    a = [9.554909596   -0.0000021389    0.0             0.0];
    e = [0.05550862    -0.000346818    -0.0000006456	0.00000000338];
    i = [2.488878       0.0025515      -0.00004903      0.000000018];
    W = [113.665524    -0.2566649      -0.00018345      0.000000357];
    P = [93.056787      0.5665496       0.00052809      0.000004882];
    mu_p = 3.79406260611373e7;
elseif Planet == 7;
    L = [314.055005      428.4669983     -0.00000486     0.000000006];
    a = [19.218446062   -0.0000000372	 0.00000000098  0];
    e = [0.04629590	 -0.000027337	 0.0000000790	 0.00000000025];
    i = [0.773196	 -0.0016869	 0.00000349	 0.000000016];
    W = [74.005947	 0.0741461	 0.00040540	 0.000000104];
    P = [173.005159	 0.0893206	 -0.00009470	 0.000000413];
    mu_p = 5.79454900707188e6;
elseif Planet == 8;
    L = [304.348665	 218.4862002	 0.00000059	 -0.000000002];
    a = [30.110386869	 -0.0000001663	 0.00000000069 0.0];
    e = [0.00898809	 0.000006408	 -0.0000000008	 -0.00000000005];
    i = [1.769952	 0.0002257	 0.00000023	 0.0];
    W = [131.784057	 -0.0061651	 -0.00000219	 -0.000000078];
    P = [48.123691	 0.0291587	 0.00007051	 -0.000000023];
    mu_p = 6.83653406387926e6;
elseif Planet == 9;
    % Pluto
    L = [238.92903833   145.20780515	0.0             0.0];
    a = [39.48211675   -0.00031596      0.0             0.0];
    e = [0.24882730     0.00005170      0.0             0.0];
    i = [17.14001206	0.00004818      0.0             0.0];
    W = [110.30393684  -0.01183482      0.0             0.0];
    P = [224.06891629  -0.04062942      0.0             0.0];
    mu_p = 9.81600887707005e2;
end

coeff.L = L;
coeff.a = a;
coeff.e = e;
coeff.i = i;
coeff.W = W;
coeff.P = P;
coeff.mu_p = mu_p;

end