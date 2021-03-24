function del_V= delV_Combined(r1,r2,r3,i1,i2,i3)
 dV =zeros(1,4);
 dV(1:2)=delV_Hohmann(r1,r2,r3);
 dV(3:4)=delV_inc(i1,i2,i3,r1,r3);
 del_V=sum(dV);
end