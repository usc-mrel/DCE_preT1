% This is the pre-contrast T1 mapping script, including conventional 
% approach (SENSE + DESPOT1) and direct estimation approach
% STARDCE package is required for DESPOT1, DCE_direct_recon (MREL) is required 
% for optimizationand Clinical_Recon_Pack (MREL) and ESP toolbox 
% (Justin Haldar and Tae Hyung Kim) are required for error evaluation.
% The script will automatically run conventional and direct T1 mapping, 
% save results and print errors for 17 undersampling levels.

clear;
clc;
close all;

%%
addpath(genpath('./SENSE+DESPOT'));
addpath(genpath('./Direct_T1'));
addpath(genpath('./minFunc_2012'));

%% Direct estimation
Direct_demo;

%% Conventional appraoch
SENSE_CG_demo;

%% Evaluation
addpath(genpath('./esp_code'));

evaluation;