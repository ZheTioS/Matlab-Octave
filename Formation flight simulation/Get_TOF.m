function [T,Orbit] = Get_TOF(a,u,e,Body_Mass,Dist_to_Sun,Norbits)
% This function will input various orbital parameters and determine the spacecraft time of flight as well as the
% conic section of the trajectory

% INPUT - a = semimajor axis (km)
        % u = gravitational parameter
        % e = eccentricity
        % Body_Mass = Mass of the orbital body (kg)
        % Dist_to_Sun = The distance from the Sun (km) at which the orbital body orbits in the heliocentric frame
        
% OUTPUT - T = Time of flight of the spacecraft (sec)
         % Orbit = Conic section of the orbit trajectory (string type)


    if e < 1 && e >= 0 % if the orbit is circular or elliptical (we know by the value of eccentricity [0.0 - 1.0])
        T = 2*pi*Norbits*sqrt(a^3/u); % sec
        if e == 0
            Orbit = 'A Circular Orbit'; % orbit is circular
        else
            Orbit = 'An Elliptical Orbit'; % orbit is elliptical
        end
    elseif e > 1 % if the orbit is hyperbolic (we know by the value of eccentricity [e>1])
        % Find the sphere of influence in relation to sun
        M_sun = 1.989e30; % Mass of the Sun in kg 
        R_soi = Dist_to_Sun*(Body_Mass/M_sun)^(2/5); % calculate the sphere of influence radius of the orbit body in km
        % find true anomaly at the edge of the sphere of influence
        R_mag = R_soi; % Distance magnitude is equal to the sphere of influence radius
        theta = acos(((a*(1-e^2)/R_mag)-1)/e); % using geometry, find the true anomaly at which the spacecraft enters the sphere of influence in radians
        
        % find hyperbolic anomaly F
        F = 2*atanh(sqrt((e-1)/(e+1))*tan(theta/2)); % using geometry, find the hyperbolic anomaly from the true anomaly in radians
        
        % find time of flight out to edge of sphere of influence
        T = sqrt((-a)^3/u)*(e*sinh(F) - F); % from periapsis to the edge of the sphere of influence
        Orbit = 'A Hyperbolic Orbit'; % orbit is hyperbolic
        
    elseif e == 1 % if the orbit is parabolic (we know by the value of eccentricity [e=1])
        % Find the sphere of influence in relation to sun
        M_sun = 1.989e30; % Mass of the Sun in kg  
        R_soi = Dist_to_Sun*(Body_Mass/M_sun)^(2/5)*a; % calculate the sphere of influence radius of the orbit body in km
        
        % find true anomaly at the sphere of influence
        R_mag = R_soi; % Distance magnitude is equal to the sphere of influence radius
        theta = acos(((a*(1-e^2)/R_mag)-1)/e); % using geometry, find the true anomaly at which the spacecraft enters the sphere of influence in radians
        


        % find parabolic anomaly D using Barker's equation
        D = tan(theta/2); % using geometry, find the parabolic anomaly from the true anomaly in radians
        % find time of flight from periapsis out to edge of sphere of influence
        T = sqrt(2*(a*(1-e^2))^3/u)*(D + D^3/3); % from periapsis to the edge of the sphere of influence
        Orbit = 'A Parabolic Orbit'; % orbit is parabolic
        
    else
        % should the there be an error in calculation and somehow the value of eccentricity is negative, throw an error message
        T = 0; % in this case there will be no time of flight
        Orbit = 'No orbit'; % no orbit
        msgbox('Negative or complex eccentricity... Reenter data!','Data Error'); % error message
    end
end