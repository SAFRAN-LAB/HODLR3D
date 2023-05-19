clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the data loaded from a given file, the code outputs the numerical
% rank for a given tolerance of 'tolr'
% OUTPUT: 'r_f' numerical rank of the face sharing interaction, 'r_e' numerical rank of the edge sharing
% interaction, 'r_v' numerical rank of the vertex sharing interaction, 'r_w' numerical rank of the well-separated interaction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = load('output_file_125_3.mat');
tolr = 1e-14;
r_f = numerical_rank(data.svd_f,tolr)
r_e = numerical_rank(data.svd_e,tolr)
r_v = numerical_rank(data.svd_v,tolr)
r_w = numerical_rank(data.svd_w,tolr)