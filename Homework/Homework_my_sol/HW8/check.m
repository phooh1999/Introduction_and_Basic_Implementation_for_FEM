clear
clc

% load('A_my.mat')
% A_my = full(A);
% load('A_check.mat')
% A_ref = full(A);
% 
% A_check = A_ref-A_my;
% 
% Mall_A = max(A_check,[],'all');

load('An_my.mat')
A_my_1 = full(A1);
A_my_2 = full(A2);
A_my_3 = full(A3);
A_my_4 = full(A4);
A_my_5 = full(A5);
A_my_6 = full(A6);
A_my_7 = full(A7);
A_my_8 = full(A8);
load('An_check.mat')
A_ref_1 = full(A1);
A_ref_2 = full(A2);
A_ref_3 = full(A3);
A_ref_4 = full(A4);
A_ref_5 = full(A5);
A_ref_6 = full(A6);
A_ref_7 = full(A7);
A_ref_8 = full(A8);

A_check_1 = A_ref_1-A_my_1;
A_check_2 = A_ref_2-A_my_2;
A_check_3 = A_ref_3-A_my_3;
A_check_4 = A_ref_4-A_my_4;
A_check_5 = A_ref_5-A_my_5;
A_check_6 = A_ref_6-A_my_6;
A_check_7 = A_ref_7-A_my_7;
A_check_8 = A_ref_8-A_my_8;

Mall_A1 = max(A_check_1,[],'all');
Mall_A2 = max(A_check_2,[],'all');
Mall_A3 = max(A_check_3,[],'all');
Mall_A4 = max(A_check_4,[],'all');
Mall_A5 = max(A_check_5,[],'all');
Mall_A6 = max(A_check_6,[],'all');
Mall_A7 = max(A_check_7,[],'all');
Mall_A8 = max(A_check_8,[],'all');


% load('b_my.mat')
% b_my = b;
% load('b_check.mat')
% b_ref = b;
% 
% b_check = b_ref-b_my;
% 
% Mall_b = max(b_check,[],'all');