function [delV_1,delV_2]=delV_Hohmann(ri,r2,rf)
	mu= 398600.4;
	delV_1=(sqrt(mu/ri)*(sqrt(2*r2/(ri+r2))-1))+(sqrt(mu/r2)*(1-sqrt(2*ri/(ri+r2))));
	delV_2=(sqrt(mu/r2)*(sqrt(2*rf/(r2+rf))-1))+(sqrt(mu/rf)*((sqrt(2*r2/(r2+rf)))-1));
	printf(’The Delta V for the first manoeuver is’);
	disp(delV_1);
	printf(’The Delta V for the second manoeuver’);
	disp(delV_2);
end