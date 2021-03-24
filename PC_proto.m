A = input('Enter data set (4x4):');


for N =1:4
    
    for B = 1:4
        pr = pr + A(N,B);
    end
    if (mod(pr,2) ~= 1)
            A(N,5)=0;
    else
            A(N,5)=1;
    end
    pr = 0;
end
for N =1:4
    pr = 0;
    for B = 1:4
        pr = pr + A(B,N);
    end
    if (mod(pr,2) ~= 1)
            A(5,N)=0;
    else (pr/2 == 0)
            A(5,N)=1;
    end
end
