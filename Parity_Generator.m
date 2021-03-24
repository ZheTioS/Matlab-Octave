function [ C ] = Parity_Generator( A )
%Parity_Generator Summary of this function goes here
%   Detailed explanation goes here
global temp
for I =1:4
    
    for J = 1:4
        temp = temp + A(I,J);
    end
    if (mod(temp,2) ~= 1)
            A(I,5)=0;
    else
            A(I,5)=1;
    end
    temp = 0;
end
for I =1:4
    temp = 0;
    for J = 1:4
        temp = temp + A(J,I);
    end
    if (mod(temp,2) ~= 1)
            A(5,I)=0;
    else (temp/2 == 0)
            A(5,I)=1;
    end
end
C = A;
end

