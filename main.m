clc; clear; close all;

addpath(genpath('src'));
addpath('config');
addpath('data');

mode = 2;

switch mode
    case 1
% 正演
main_forward();
    case 2
% 反演
main_inversion();
end