%F is the Failure Producing Value
%A is the Tensile Strength\
%r is our calculated radius
%rvec is a range of possible radii
%Smin - Minimum Tensile Strength
%Smax - Maxium Tensile Strength
%kmin - safety min
%kmax - safety max
%kvecmax - max range of safety values
%kvecmin - min range of safety values

DHE = .002489;
DMAT = 920;

%r = 24.876; %m^3
rvec = linspace(20,30,10);
t = 2.032*10^-5; %thickness in meters
[T, a, P, rho] = atmoscoesa(30000);
Smin = 19.31*10^6; %Pascals / Meter
Smax = 22.06*10^6; %Pascals / Meter

k = linspace(1,2,10);

r = (1000./(rho*(4/3)*pi - (DMAT*4*10*pi.*k)./(2*Smax) - DHE*(4/3)*pi)).^(1/3); %Original Material 105F
r2 = (1000./(rho*(4/3)*pi - (918*4*10*pi.*k)./(2*1200*6894.76) - DHE*(4/3)*pi)).^(1/3); %1052
r3 = (1000./(rho*(4/3)*pi - (925*4*10*pi.*k)./(2*6894.76*4300) - DHE*(4/3)*pi)).^(1/3); %108A
figure
plot(k,r,k,r2,k,r3)
legend('Polyetheylene Bapolene 105F','Polyetheylene Bapolene 1052','Polyetheylene Bapolene 108A');
title('Safety Factor vs Radius')
xlabel('Safety Factor')
ylabel('Radius in meters')

 fopen('results.txt','w')
 fprintf(fid,'-------------------')
 fprintf(fid,'\r\n')
 fprintf(fid,'Polyetheylene Bapolene 105F Radius Values');
 fprintf(fid,'\r\n')
 fprintf(fid,'-------------------')
 fprintf(fid,'\r\n')
 for i = 1:length(r)
     fprintf(fid,'%.2f',r(i)) 
     fprintf(fid,'\r\n');
 end
 fprintf(fid,'-------------------')
 fprintf(fid,'\r\n');
 fprintf(fid,'Polyetheylene Bapolene 108A Radius Values');
 fprintf(fid,'\r\n')
 fprintf(fid,'-------------------')
 fprintf(fid,'\r\n')
  for i = 1:length(r2)
     fprintf(fid,'%.2f',r3(i)) 
     fprintf(fid,'\r\n');
 end
 fprintf(fid,'-------------------')
 fprintf(fid,'\r\n');
 fprintf(fid,'Polyetheylene Bapolene 105F Radius Values');
 fprintf(fid,'\r\n')
 fprintf(fid,'-------------------')
 fprintf(fid,'\r\n')
  for i = 1:length(r3)
     fprintf(fid,'%.2f',r2(i)) 
     fprintf(fid,'\r\n');
 end
 fprintf(fid,'-------------------')
fprintf(fid,'\r\n');

