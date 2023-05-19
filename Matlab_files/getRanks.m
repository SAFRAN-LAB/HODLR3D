clc;
clear all;
close all;
data = load('output_file_125_3.mat');
tolr = 1e-14;
r_f = numerical_rank(data.svd_f,tolr)
r_e = numerical_rank(data.svd_e,tolr)
r_v = numerical_rank(data.svd_v,tolr)
r_w = numerical_rank(data.svd_w,tolr)