function [ sum, invA, invB, traA, traB, product ] = matops( A, B )
%matops Takes 2 matrices A and B and returns sum,product,inverse and
%transpose of the matrices.
%   sum - Sum of A and B
%   invA - inverse of A
%   invB - inverse of B
%   traA - transpose of A
%   traB - transpose of B
%   product - Product of A abd B
sum = A+B
invA= inv(A)
invB = inv(B)
traA = A'
traB = B'
product = A*B

end

