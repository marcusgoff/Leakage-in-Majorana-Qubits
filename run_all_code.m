%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: run_all_code.m
% Purpose: Generates all figures from our paper 
%          "Leakage at Zero Temperature from Changes in Chemical Potential
%           in Majorana Qubits"
%
% Note: you may alternatevly generate individual figures by opening and 
% running their corresponding ".m" files.
%
% Author: Marcus C. Goffage
% Affiliation: University of New South Wales
%
% Paper: "Decoherence in Majorana Qubits by 1/f Noise"
% Paper Authors: M. C. Goffage^1, A. Alase^2, M. C. Cassidy^1, 
%                S. N. Coppersmith^1
% Affiliations: ^1 University of New South Wales
%               ^2 University of Sydney 
%               
%
%----------------------------------------%

run run_fig_2_main.m
run run_fig_2_inset.m

run 'Supplement Plots'/run_fig_S2.m
run 'Supplement Plots'/run_fig_S3.m
run 'Supplement Plots'/run_fig_S4.m
run 'Supplement Plots'/run_fig_S5.m