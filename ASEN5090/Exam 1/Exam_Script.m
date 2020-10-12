clear all;close all;clc
c   = 2.99792458e8;    % GPS acceptd speed of light, m/s
%Create equator and near north pole user positions
EQUAECEF = 1000*[6380 0 0];
POLEECEF = 1000*[0 0 6380];
userpos = [EQUAECEF; POLEECEF];

%Create stand-in almanac entry to be used for exam
[gps_ephem,gps_ephem_cell] = read_GPSyuma('YUMA245.ALM',2);
ephem = gps_ephem(1,:);
ephem(4) = 0; %set eccen. to 0 for circ orbit
ephem(5) = sqrt(26560000); %Set square root of semi major axis as square root of circ. orbit radius
ephem(6) = 0;
ephem(7) = deg2rad(55); %Set inc, as exactly 55 degress
ephem(9) = 0; %Change rate of change or right ascension to 0 to account for assumption Eearth is not rotating.

% tvec = gps_ephem_cell{1,1}.Toe:60:gps_ephem_cell{1,1}.Toe+24*3600;
tvec = 0:30:86400;
tvec = tvec + 2*86400;

[health,pos] = broadcast2pos_alt(ephem,[2121*ones(length(tvec),1) tvec'],1);

for j = 1:2
    for i = 1:size(pos,1)
        [az(i,j),el(i,j),range(i,j)] = compute_azelrange_alt(userpos(j,:),pos(i,:));
    end
end
    figure
    for j = 1:size(pos,1)
        if el(j,1)>0
            plotAzEl(az(j,1),el(j,1),0)
        end
    end
    title('Equator Observer Sky Plot')
    figure
    for j = 1:size(pos,1)
        if el(j,2)>0
            plotAzEl(az(j,2),el(j,2),0)
        end
    end
    title('North Pole Observer Sky Plot')
%%
%Determine min range vals.
EQ = 26560-6380;
temp = max(el(:,2));
temp = find(el(:,2) == temp);
Pole = range(temp,2);
%range at beginning of pass found when elevation first goes positive
EQ_passrange = range(735,1);
Pole_passrange = range(1107,2);
%TOF of range computed in part c
EQ_time = EQ_passrange*(1/c);
Pole_time = Pole_passrange*(1/c);
%max elevation for both observer positions
maxel = max(el);
%approximate duration of satellite pass
EQ_passtime = (tvec(1341) - tvec(735))/3600;
Pole_passtime = (tvec(1688) - tvec(1107))/3600;
%calculating max and min doppler shifts for obeserver positions
%minimum doppler for L1 for both positions is 0 because the range rate will
%always change from positive to negative at some point for both positions,
%nad this indicates a change from positive to negative doppler shift, so
%the minimum must be 0
ft = 1575.42;
rdotEQ = diff(range(:,1));
rdotPole = diff(range(:,2));
maxfdEQ = (-max(rdotEQ)/c)*ft;
maxfdPole = (-max(rdotPole)/c)*ft;
figure
for i = 1:numel(tvec)
    if el(i,1)
    plot(tvec-2*86400,el(:,1))
    end
end