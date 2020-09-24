%Body temp plot with error
x = 1:10;
body_temp_volt = [1340 1367 1389 1396 1401 1408 1330 1387 1401 1400]; %Copy and paste your 10 body temp voltage measurements here
body_temp = [];

for i = 1:length(body_temp_volt)
    body_temp = [body_temp voltToTemp(body_temp_volt(i))];
end

error_volt = 2;%Use whatever observed error you had for your 10 body temp measurements
error = [];
for i = 1:length(body_temp_volt)
    error = [error voltToTemp(error_volt)];
end
%The following weighted mean calculation is done assuming weight is equal for all ten measurements
A = body_temp / 10;
Weighted_mean_body_temp = sum(A);
Weighted_mean_body_temp = repelem(Weighted_mean_body_temp,length(x));
%Standard deviation calculation
body_temp_mean = (sum(body_temp)/10);
B = (body_temp - body_temp_mean).^2;
C = (sum(B)/10);
D = sqrt(C);
error_of_mean = (D/sqrt(length(x)));
error_of_mean = repelem(error_of_mean,length(x));
hold on

errorbar(x,body_temp,error);
plot(x,Weighted_mean_body_temp);

plot(x,(Weighted_mean_body_temp + error_of_mean),'--g');
plot(x,Weighted_mean_body_temp - error_of_mean,'--r');
xlabel('Measurement Number');
ylabel('Degrees Celsius');
title('Body Temperature in Degrees Celsius-Boylston/Herrin');


legend('Body temp in degrees C','Weighted Mean Body Temperature','Weighted Mean Body Temp + Error','Weighted Mean Body Temp - Error');

