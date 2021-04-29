%This file serves as the main file should the grader wish run through the
%functions and scripts used to complete this project. This main script
%serves to run each of the subsequent main scripts used throughout the
%project, and will make a best attempt at making it easy for the grader to
%reproduce results and plots shown in the report.
%% Problem 1 -- Dynamics and Measurement Verification
Test
%% Problem 2 -- Known Model and Truth
CKF_SNC
%% Problem 3 -- Unknown Issues and no Truth
CKF_SNC_p3
%{
  As a note, several of the plots that are presented in the report were 
  generated using multiple runs of one of these filter in 
  conjunction with saved workspaces that were all loaded and used in one
  instance to give the desired results, e.g. the plot of 4 B-plane
  estimates and their 3-sigma covariance ellipses
%}