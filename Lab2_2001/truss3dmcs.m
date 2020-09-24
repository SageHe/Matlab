    function truss3dmcs(inputfile)
% function truss2d(inputfile,outputfile)
%
% Stochastic analysis of 2-D statically determinate truss by
% Monte Carlo Simulation. Only positions and strength of joints 
% treated as random variables
%
% Assumption: variation of joint strength and positions described 
%             via Gaussian distributions
% 
%             joint strength : mean = 15
%                              coefficient of varation = 0.1
%             joint position : 
%                              coefficient of varation = 0.01
%                              (defined wrt to maximum dimension of truss)
%
%             number of samples is set to 1e4
%
% Input:  inputfile  - name of input file
%
% Author: Kurt Maute for ASEN 2001, Oct 13 2012

% parameters
jstrmean   = 4.8;    % mean of joint strength
jstrcov    = .4;   % coefficient of variation of joint strength
jposcov    = .01;  % coefficient of variation of joint position
numsamples = 1e5;   % number of samples

% read input file
[joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs]=readinput(inputfile);

% determine extension of truss
ext_x=max(joints(:,1))-min(joints(:,1));   % extension in x-direction
ext_y=max(joints(:,2))-min(joints(:,2));   % extension in y-direction
ext_z=max(joints(:,3))-min(joints(:,3));   % extension in z-direction
ext  =max([ext_x,ext_y,ext_z]);

% loop overall samples
numjoints=size(joints,1);       % number of joints
maxforces=zeros(numsamples,1);  % maximum bar forces for all samples
maxreact=zeros(numsamples,1);   % maximum support reactions for all samples
failure=zeros(numsamples,1);    % failure of truss

for is=1:numsamples 
    
    % generate random joint strength limit
    varstrength = (jstrcov*jstrmean)*randn(1,1);
    
    jstrength = jstrmean + varstrength;
    
    % generate random samples
    varjoints = (jposcov*ext)*randn(numjoints,3); %originally (numjoints,2)
    
    % perturb joint positions
    randjoints = joints + varjoints;
    
    % compute forces in bars and reactions
    [barforces,reacforces] = forceanalysis(randjoints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs);
    
    % determine maximum force magnitude in bars and supports
    maxforces(is) = max(abs(barforces));
    maxreact(is)  = max(abs(reacforces));
    
    % determine whether truss failed
    failure(is) = maxforces(is) > jstrength || maxreact(is) > jstrength;
end

figure(1);
subplot(1,2,1);
hist(maxforces,30);
title('Histogram of maximum bar forces');
xlabel('Magnitude of bar forces');
ylabel('Frequency');

subplot(1,2,2);
hist(maxreact,30);
title('Histogram of maximum support reactions');
xlabel('Magnitude of reaction forces');
ylabel('Frequency');

fprintf('\nFailure probability : %e \n\n',sum(failure)/numsamples);

end

function writeoutput(outputfile,inputfile,barforces,reacforces,joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs)
% writeoutput(outputfile,inputfile,barforces,reacforces,joints,connectivity,reacjoints,reacvecs,loadjoints,loadvecs);
%
% output analysis results
%
% input:  outputfile   - name of output file
%         inputfile    - name of input file
%         barforces    - force magnitude in bars
%         reacforces   - reaction forces
%         joints       - coordinates of joints
%         connectivity - connectivity 
%         reacjoints   - joint id where reaction acts on
%         reacvecs     - unit vector associated with reaction force
%         loadjoints   - joint id where external load acts on
%         loadvecs     - load vector
%
%
% Author: Kurt Maute, Sept 21 2011

% open output file
fid=fopen(outputfile,'w');

% write header
fprintf(fid,'2-D Truss analysis\n');
fprintf(fid,'------------------\n\n');
fprintf(fid,'Date: %s\n\n',datestr(now));

% write name of input file
fprintf(fid,'Input file: %s\n\n',inputfile);

% write coordinates of joints
fprintf(fid,'Joints:         Joint-id  x-coordinate y-coordinate\n');
for i=1:size(joints,1)
    fprintf(fid,'%17d %12.2f %12.2f\n',i,joints(i,1),joints(i,2));
end
fprintf(fid,'\n\n');

% write external loads
fprintf(fid,'External loads: Joint-id  Force-x      Force-y\n');
for i=1:size(loadjoints,1)
    fprintf(fid,'%17d %12.2f %12.2f\n',loadjoints(i),loadvecs(i,1),loadvecs(i,2));
end
fprintf(fid,'\n');
    
% write connectivity and forces
fprintf(fid,'Bars:           Bar-id    Joint-i      Joint-j     Force    (T,C)\n');
for i=1:size(connectivity,1)
    if barforces(i)>0;tc='T';else tc='C';end
    fprintf(fid,'%17d   %7d %12d    %12.3f     (%s)\n',...
        i,connectivity(i,1),connectivity(i,2),abs(barforces(i)),tc);
end
fprintf(fid,'\n');

% write connectivity and forces
fprintf(fid,'Reactions:      Joint-id  Uvec-x       Uvec-y      Force\n');
for i=1:size(reacjoints,1)
    fprintf(fid,'%17d %12.2f %12.2f %12.3f\n',reacjoints(i),reacvecs(i,1),reacvecs(i,2),reacforces(i));
end

% close output file
fclose(fid);

end