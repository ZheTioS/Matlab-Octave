function [delV_3,delv_4]=delV_inc(i1,i2,i3,ri,rf)
 mu=398600.4;
 del_i1=i1-i2;
 del_i2=(i2-i3)*(-1);
 delV_3=2*sqrt(mu/ri)*sind(del_i1/2);
 delV_4=2*sqrt(mu/rf)*sind(del_i2/2);
 printf(’The Delta V for the first inclination change is’);
 disp(delV_3);
 printf(’The Delta V for the second inclination change is’);
 disp(delV_4);
end