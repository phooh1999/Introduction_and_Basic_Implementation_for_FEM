clear
clc

% load('A_my.mat')
% A_my = full(A);
% load('A_ref.mat')
% A_ref = full(A);
% 
% A_check = A_ref-A_my;
% 
% Mall_A = max(A_check,[],'all');

load('b_my.mat')
b_my = b;
load('b_check.mat')
b_ref = b;

b_check = b_ref-b_my;

Mall_b = max(b_check,[],'all');