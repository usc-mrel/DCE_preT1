% This is the pre-contrast T1 mapping script, including conventional 
% approach (SENSE + DESPOT1) and direct estimation approach
% STARDCE package is required for DESPOT1, DCE_direct_recon (MREL) is required 
% for optimizationand Clinical_Recon_Pack (MREL) and ESP toolbox 
% (Justin Haldar and Tae Hyung Kim) are required for error evaluation.
% The script will automatically run conventional and direct T1 mapping, 
% save results and print errors

clear;
clc;
close all;

%%
addpath(genpath('./SENSE+DESPOT'));
addpath(genpath('./Direct_T1'));
addpath(genpath('./STARDCE-master-NEW'));
addpath(genpath('./DCE_direct_recon'));
addpath(genpath('./esp_code'));

%% Conventional appraoch
SENSE_CG;

%% Direct estimation
Direct;

%% Error evaluation
addpath(genpath('./Clinical_Recon_Pack-master'));
evaluation;