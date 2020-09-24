C_L = ones(1,32);
C_D = ones(1,32);
Data = xlsread('AirfoilPressure_S011_G03.xlsx');
c = 0.0889;
%Average data before C_p calc1
for x = 1:28
    dummy = Data(:,x);
    for i = 1:12
        airspeedmeans(i,1) = mean(dummy((i-1)*20 + 1 : i*20,:));
    end
    airspeedmeans(1:12,x) = airspeedmeans(:,1);
end
valvep = airspeedmeans(:,7:22);

scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];

scanivalve_pos = scanivalve_pos*0.0254;

[Cn, Ca, CP] = FindCp(airspeedmeans,scanivalve_pos);

[CL, CD] = draglift(airspeedmeans,Ca,Cn);

CL = CL*-1;
CD = CD*-1;

scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 3.5 0; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];
scanivalve_pos = scanivalve_pos*0.0254;
norm_chord_length = (scanivalve_pos(:,1)/c);

hold on

for i = 1:12
    figure(1)
    plot(norm_chord_length,CP(i,:))
%     hold on
    axis ij
    title('Pressure Coefficient VS Normalized Chord Length-Group 3')
    xlabel('Normaliz Chord Length (X/C)')
    ylabel('Pressure Coefficient Cp')
end
    legend('9 M/S, -8 degrees','17 M/S, -8 degrees','34M/S -8 degrees','9 M/S, 0 degrees','17 M/S, 0 degrees','34M/S 0 degrees','9 M/S, 0 degrees','17 M/S, 0 degrees','34M/S 0 degrees','9 M/S, 16 degrees','17 M/S, 16 degrees','34M/S 16 degrees')
    

% hold off

attack_angle = airspeedmeans([1 4 7 10],23);

figure(3)
plot(attack_angle,CD([2 5 8 11]))
hold on

figure(4)
plot(attack_angle,CL([2 5 8 11]))
hold on
% function [CL CD] = draglift(airspeedmeans,Ca,Cn)
% for x = 1:12
%     CL(1,x) = Cn(x)*cosd(airspeedmeans(x,23)) - Ca(x)*sind(airspeedmeans(x,23));
%     CD(1,x) = Cn(x)*sind(airspeedmeans(x,23)) + Ca(x)*cosd(airspeedmeans(x,23));
% end
% end
C_L([attack_angle + 16]) = CL([2 5 8 11]);
C_D([attack_angle + 16]) = CD([2 5 8 11]);

Data = xlsread('AirfoilPressure_S011_G09.xlsx');
c = 0.0889;
%Average data before C_p calc1
for x = 1:28
    dummy = Data(:,x);
    for i = 1:12
        airspeedmeans(i,1) = mean(dummy((i-1)*20 + 1 : i*20,:));
    end
    airspeedmeans(1:12,x) = airspeedmeans(:,1);
end
valvep = airspeedmeans(:,7:22);

scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];

scanivalve_pos = scanivalve_pos*0.0254;

[Cn, Ca, CP] = FindCp(airspeedmeans,scanivalve_pos);

[CL, CD] = draglift(airspeedmeans,Ca,Cn);

CL = CL*-1;
CD = CD*-1;

scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 3.5 0; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];
scanivalve_pos = scanivalve_pos*0.0254;
norm_chord_length = (scanivalve_pos(:,1)/c);

hold on

for i = 1:12
    figure(2)
    plot(norm_chord_length,CP(i,:))
    hold on
    axis ij
    title('Pressure Coefficient VS Normalized Chord Length-Group 9')
    xlabel('Normaliz Chord Length (X/C)')
    ylabel('Pressure Coefficient Cp')
end

    legend('9 M/S, -12 degrees','17 M/S, -12 degrees','34M/S -12 degrees','9 M/S, -4 degrees','17 M/S, -4 degrees','34M/S -4 degrees','9 M/S, 4 degree','17 M/S, 4 degree','34M/S 4 degree','9 M/S, 12 degrees','17 M/S, 12 degrees','34M/S 12 degrees')

attack_angle = airspeedmeans([1 4 7 10],23);

figure(3)
plot(attack_angle,CD([2 5 8 11]))
hold on

figure(4)
plot(attack_angle,CL([2 5 8 11]))
hold on

C_L([attack_angle + 16]) = CL([2 5 8 11]);
C_D([attack_angle + 16]) = CD([2 5 8 11]);

Data = xlsread('AirfoilPressure_S011_G05.xlsx');
c = 0.0889;
%Average data before C_p calc1
for x = 1:28
    dummy = Data(:,x);
    for i = 1:12
        airspeedmeans(i,1) = mean(dummy((i-1)*20 + 1 : i*20,:));
    end
    airspeedmeans(1:12,x) = airspeedmeans(:,1);
end
valvep = airspeedmeans(:,7:22);

scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];

scanivalve_pos = scanivalve_pos*0.0254;

[Cn, Ca, CP] = FindCp(airspeedmeans,scanivalve_pos);

[CL, CD] = draglift(airspeedmeans,Ca,Cn);

CL = CL*-1;
CD = CD*-1;

scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 3.5 0; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];
scanivalve_pos = scanivalve_pos*0.0254;
norm_chord_length = (scanivalve_pos(:,1)/c);

hold on

for i = 1:12
    figure(1)
    plot(norm_chord_length,CP(i,:))
    hold on
    axis ij
    title('Pressure Coefficient VS Normalized Chord Length')
    xlabel('Normaliz Chord Length (X/C)')
    ylabel('Pressure Coefficient Cp')
end

% hold off

attack_angle = airspeedmeans([1 4 7 10],23);

figure(3)
plot(attack_angle,CD([2 5 8 11]))
hold on

figure(4)
plot(attack_angle,CL([2 5 8 11]))
hold on

C_L([attack_angle + 16]) = CL([2 5 8 11]);
C_D([attack_angle + 16]) = CD([2 5 8 11]);

Data = xlsread('AirfoilPressure_S011_G07.xlsx');
c = 0.0889;
%Average data before C_p calc1
for x = 1:28
    dummy = Data(:,x);
    for i = 1:12
        airspeedmeans(i,1) = mean(dummy((i-1)*20 + 1 : i*20,:));
    end
    airspeedmeans(1:12,x) = airspeedmeans(:,1);
end
valvep = airspeedmeans(:,7:22);

scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];

scanivalve_pos = scanivalve_pos*0.0254;

[Cn, Ca, CP] = FindCp(airspeedmeans,scanivalve_pos);

[CL, CD] = draglift(airspeedmeans,Ca,Cn);

CL = CL*-1;
CD = CD*-1;

scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 3.5 0; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];
scanivalve_pos = scanivalve_pos*0.0254;
norm_chord_length = (scanivalve_pos(:,1)/c);

hold on

for i = 1:12
    figure(1)
    plot(norm_chord_length,CP(i,:))
    hold on
    axis ij
    title('Pressure Coefficient VS Normalized Chord Length')
    xlabel('Normaliz Chord Length (X/C)')
    ylabel('Pressure Coefficient Cp')
end

% hold off

attack_angle = airspeedmeans([1 4 7 10],23);

figure(3)
plot(attack_angle,CD([2 5 8 11]))
hold on

figure(4)
plot(attack_angle,CL([2 5 8 11]))
hold on

C_L([attack_angle + 16]) = CL([2 5 8 11]);
C_D([attack_angle + 16]) = CD([2 5 8 11]);

Data = xlsread('AirfoilPressure_S011_G01.xlsx');
c = 0.0889;
%Average data before C_p calc1
for x = 1:28
    dummy = Data(:,x);
    for i = 1:12
        airspeedmeans(i,1) = mean(dummy((i-1)*20 + 1 : i*20,:));
    end
    airspeedmeans(1:12,x) = airspeedmeans(:,1);
end
valvep = airspeedmeans(:,7:22);

scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];

scanivalve_pos = scanivalve_pos*0.0254;

[Cn, Ca, CP] = FindCp(airspeedmeans,scanivalve_pos);

[CL, CD] = draglift(airspeedmeans,Ca,Cn);

CL = CL*-1;
CD = CD*-1;

scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 3.5 0; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];
scanivalve_pos = scanivalve_pos*0.0254;
norm_chord_length = (scanivalve_pos(:,1)/c);

hold on

for i = 1:12
    figure(1)
    plot(norm_chord_length,CP(i,:))
    hold on
    axis ij
    title('Pressure Coefficient VS Normalized Chord Length')
    xlabel('Normaliz Chord Length (X/C)')
    ylabel('Pressure Coefficient Cp')
end

% hold off

attack_angle = airspeedmeans([1 4 7 10],23);

figure(3)
plot(attack_angle,CD([2 5 8 11]))
hold on

figure(4)
plot(attack_angle,CL([2 5 8 11]))
hold on

C_L([attack_angle + 16]) = CL([2 5 8 11]);
C_D([attack_angle + 16]) = CD([2 5 8 11]);

Data = xlsread('AirfoilPressure_S011_G11.xlsx');
c = 0.0889;
%Average data before C_p calc1
for x = 1:28
    dummy = Data(:,x);
    for i = 1:12
        airspeedmeans(i,1) = mean(dummy((i-1)*20 + 1 : i*20,:));
    end
    airspeedmeans(1:12,x) = airspeedmeans(:,1);
end
valvep = airspeedmeans(:,7:22);

scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];

scanivalve_pos = scanivalve_pos*0.0254;

[Cn, Ca, CP] = FindCp(airspeedmeans,scanivalve_pos);

[CL, CD] = draglift(airspeedmeans,Ca,Cn);

CL = CL*-1;
CD = CD*-1;

scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 3.5 0; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];
scanivalve_pos = scanivalve_pos*0.0254;
norm_chord_length = (scanivalve_pos(:,1)/c);

hold on

for i = 1:12
    figure(1)
    plot(norm_chord_length,CP(i,:))
    hold on
    axis ij
    title('Pressure Coefficient VS Normalized Chord Length')
    xlabel('Normaliz Chord Length (X/C)')
    ylabel('Pressure Coefficient Cp')
end

% hold off

attack_angle = airspeedmeans([1 4 7 10],23);

figure(3)
plot(attack_angle,CD([2 5 8 11]))
hold on

figure(4)
plot(attack_angle,CL([2 5 8 11]))
hold on

C_L([attack_angle + 16]) = CL([2 5 8 11]);
C_D([attack_angle + 16]) = CD([2 5 8 11]);

Data = xlsread('AirfoilPressure_S011_G13.xlsx');
c = 0.0889;
%Average data before C_p calc1
for x = 1:28
    dummy = Data(:,x);
    for i = 1:12
        airspeedmeans(i,1) = mean(dummy((i-1)*20 + 1 : i*20,:));
    end
    airspeedmeans(1:12,x) = airspeedmeans(:,1);
end
valvep = airspeedmeans(:,7:22);

scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];

scanivalve_pos = scanivalve_pos*0.0254;

[Cn, Ca, CP] = FindCp(airspeedmeans,scanivalve_pos);

[CL, CD] = draglift(airspeedmeans,Ca,Cn);

CL = CL*-1;
CD = CD*-1;

scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 3.5 0; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];
scanivalve_pos = scanivalve_pos*0.0254;
norm_chord_length = (scanivalve_pos(:,1)/c);

hold on

for i = 1:12
    figure(1)
    plot(norm_chord_length,CP(i,:))
    hold on
    axis ij
    title('Pressure Coefficient VS Normalized Chord Length')
    xlabel('Normaliz Chord Length (X/C)')
    ylabel('Pressure Coefficient Cp')
end

% hold off

attack_angle = airspeedmeans([1 4 7 10],23);

figure(3)
plot(attack_angle,CD([2 5 8 11]))
hold on

figure(4)
plot(attack_angle,CL([2 5 8 11]))
hold on

C_L([attack_angle + 16]) = CL([2 5 8 11]);
C_D([attack_angle + 16]) = CD([2 5 8 11]);

Data = xlsread('AirfoilPressure_S011_G15.xlsx');
c = 0.0889;
%Average data before C_p calc1
for x = 1:28
    dummy = Data(:,x);
    for i = 1:12
        airspeedmeans(i,1) = mean(dummy((i-1)*20 + 1 : i*20,:));
    end
    airspeedmeans(1:12,x) = airspeedmeans(:,1);
end
valvep = airspeedmeans(:,7:22);

scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];

scanivalve_pos = scanivalve_pos*0.0254;

[Cn, Ca, CP] = FindCp(airspeedmeans,scanivalve_pos);

[CL, CD] = draglift(airspeedmeans,Ca,Cn);

CL = CL*-1;
CD = CD*-1;

scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 3.5 0; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];
scanivalve_pos = scanivalve_pos*0.0254;
norm_chord_length = (scanivalve_pos(:,1)/c);

hold on

for i = 1:12
    figure(1)
    plot(norm_chord_length,CP(i,:))
    hold on
    axis ij
    title('Pressure Coefficient VS Normalized Chord Length')
    xlabel('Normaliz Chord Length (X/C)')
    ylabel('Pressure Coefficient Cp')
end

% hold off

attack_angle = airspeedmeans([1 4 7 10],23);

figure(3)
plot(attack_angle,CD([2 5 8 11]))
hold on

figure(4)
plot(attack_angle,CL([2 5 8 11]))
hold on

C_L([attack_angle + 16]) = CL([2 5 8 11]);
C_D([attack_angle + 16]) = CD([2 5 8 11]);


function [CL CD] = draglift(airspeedmeans,Ca,Cn)
for x = 1:12
    CL(1,x) = Cn(x)*cosd(airspeedmeans(x,23)) - Ca(x)*sind(airspeedmeans(x,23));
    CD(1,x) = Cn(x)*sind(airspeedmeans(x,23)) + Ca(x)*cosd(airspeedmeans(x,23));
end
end



