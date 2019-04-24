function [  ] = plotting( x,y,r )
%Plotting A function to plot a circle,a sphere, and y=x line
%   If coencentric 
c = input('Enter 1 - Circle \n 2- Sphere \n 3- y=x plot');
if( c == 1)
    T = 0:0.01:2*pi;
    xT = r*cos(T) + x;
    yT = r*sin(T) + y;
    plot(xT,yT)
    hold on
elseif( c == 2)
        sphere(r)
elseif ( c ==
end

