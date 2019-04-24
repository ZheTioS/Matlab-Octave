function [ matcat ] = matconcat( A,B,C )
%matconcat A function to concatnate 3 matrices A,B,C to create 
% a 3x3x3 matrix

%% Matrices A,B,C must be 3x3 2-D matrices

c = {'numeric'}; % Defining class
att = {'size',[3,3]}; % Defining dimension
validateattributes(A,c,att)
validateattributes(B,c,att)
validateattributes(C,c,att)
%% Concatenation of Matrices to form required 3X3X3 matrix 

matcat = cat(3,A,B,C);


end

