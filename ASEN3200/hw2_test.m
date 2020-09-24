%Question 1
info = [10 1 1 1;10 -1 -1 -1;8 4 -4 4;8 -2 2 -2;12 3 -3 -3;12 -3 3 3];
I = zeros(3,3);
for i = 1:size(info,1)
    I_temp = zeros(3,3);
    I_temp = [info(i,1)*(info(i,3)^2 + info(i,4)^2) -info(i,1)*info(i,2)*info(i,3) -info(i,1)*info(i,2)*info(i,4);...
         -info(i,1)*info(i,2)*info(i,3) info(i,1)*(info(i,2)^2 + info(i,4)^2) -info(i,1)*info(i,3)*info(i,4);...
         -info(i,1)*info(i,2)*info(i,4) -info(i,1)*info(i,3)*info(i,4) info(i,1)*(info(i,2)^2 + info(i,3)^2)];
     
    I = I + I_temp;
end

I_central = [60*((-4/15)^2 + (4/15)^2) -60*(4/15)*(-4/15) -60*(4/15)*(4/15);...
            -60*(4/15)*(-4/15) 60*((4/15)^2 + (4/15)^2) -60*(-4/15)*(4/15);...
            -60*(4/15)*(4/15) -60*(-4/15)*(4/15) 60*((4/15)^2 + (-4/15)^2)];

I = I - I_central;
%Question 2&3
[V,D] = eig(I);
%Question 4
R_tild = [0 0 4;0 0 0;-4 0 0];
I_cube = eye(3,3)*2;
I_sphere = eye(3,3)*(4/5);
I = I_cube + I_sphere - (6/5)*R_tild*R_tild;
%Question 6
I_one = [((.15)*(.1^2))/12 0 0; 0 ((.15)*(.1^2))/12 0;0 0 ((.15)*(.1^2))/6];
I_two = I_one;
r_tild1 = [0 -.1 0;.1 0 0;0 0 0];

I_onetwo = I_one + I_two - ((.15^2)/(.3))*r_tild1*r_tild1;

I_three = [((.15)*(.1^2))/12 0 0;0 ((.15)*(.1^2))/6 0; 0 0 ((.15)*(.1^2))/12];
I_four = I_three;

r_tild2 = [0 0 -.1; 0 0 0;.1 0 0];

I_threefour = I_three + I_four - ((.15^2)/(.3))*r_tild2*r_tild2;

I_five = [((.15)*(.1^2))/6 0 0;0 ((.15)*(.1^2))/12 0; 0 0 ((.15)*(.1^2))/12];
I_six = I_five;
r_tilda3 = [0 0 0; 0 0 -.1;0 .1 0];

I_fivesix = I_five + I_six - ((.15^2)/.3)*r_tilda3*r_tilda3;
