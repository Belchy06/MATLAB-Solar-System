function [p, v] = solarsystem(p, v, mass, stop_time, hide_animation)
% SOLARSYSTEM Solve the gravitational dynamics of n bodies.
%
% SOLARSYSTEM(p, v, mass, stop_time) receives the initial position, initial
% velocity, and mass, and simulate for stop_time seconds. The inputs p and
% v are N-by-2 matrices where the Nth row refers to the Nth celestial
% body. The two columns are the X and Y values of position and velocity,
% respectively. The input mass is an N-by-1 vector where the Nth item is
% the mass of the Nth celestial body.
%
% [p, v] = SOLARSYSTEM(...) returns the final position and final velocity
% of each object in the same format as the input. This will be used during
% marking to test the accuracy of your program.
%
% SOLARSYSTEM(..., hide_animation) is an optional flag to indicate whether
% the graphical animation may be hidden. This is an advanced feature for
% students who are aiming for a high level of achievement. During marking,
% when the computation speed of your program is being tested, it will be
% run with the hide_animation flag set to true. This allows your program to
% skip the drawing steps and run more quickly. 

    
if nargin < 5
    hide_animation = false;
end

% Physical constants
% Time constants
delta_t = 3600;

% Universal Gravitational Constant
G = 6.673 * 10^-11;
    
% Graph Variables
pens = [animatedline('LineWidth', 2),animatedline('LineWidth', 2),animatedline('LineWidth', 2),animatedline('LineWidth', 2),animatedline('LineWidth', 2)];
set(gca,'XLim', [-183610000000 183610000000],'YLim', [-183610000000 183610000000],'ZLim', [-183610000000 183610000000]);  
    
% Write your code here
% Loop through code for each second
for t=1:stop_time
    
    for planetA=1:length(p)
    
        % Retrieve data from supplied data
        planetAPosition = [p(planetA,1), p(planetA,2), p(planetA,3)];
        planetAVelocity = [v(planetA,1), v(planetA,2), v(planetA,3)];
        planetAMass = mass(planetA);

        for planetB=(planetA+1):length(p)
        
            % Retrieve data from supplied data
            planetBPosition = [p(planetB,1), p(planetB,2), p(planetB,3)];
            planetBMass = mass(planetB);

            % Calculate necessary values as defined in assignment sheet
            distanceVector = planetBPosition - planetAPosition;
            distanceVectorMagnitude = sqrt(distanceVector(1)^2 + distanceVector(2)^2 + distanceVector(3)^2);

            gravityVector = ((G * planetAMass * planetBMass) / abs(distanceVectorMagnitude^3)) * distanceVector;
            
            % Gravity Vector from planet B -> A
            % Positive value as specified in assginment sheet
            gravityVectorArray.x(planetB,planetA) = gravityVector(1);
            gravityVectorArray.y(planetB,planetA) = gravityVector(2);
            gravityVectorArray.z(planetB,planetA) = gravityVector(3);

            % Gravity Vector from planet A -> B
            % Negative value can be achieved by multiplying gravity vector component by -1
            % as the direction is just reversed. No need to recalculate components
            gravityVectorArray.x(planetA,planetB) = gravityVector(1) * -1;
            gravityVectorArray.y(planetA,planetB) = gravityVector(2) * -1;
            gravityVectorArray.z(planetA,planetB) = gravityVector(3) * -1;            

        end

        % Total Gravity in xyz components that affects planetA
        totalGravity.x = sum(gravityVectorArray.x(:,planetA));
        totalGravity.y = sum(gravityVectorArray.y(:,planetA));
        totalGravity.z = sum(gravityVectorArray.z(:,planetA));
        
        % Move total gravity into an array of acceleration values
        gravitationalAcceleration = [totalGravity.x/planetAMass totalGravity.y/planetAMass totalGravity.z/planetAMass];
        newPlanetAVelocity = planetAVelocity + (t * gravitationalAcceleration);
        newPlanetAPosition = planetAPosition + (newPlanetAVelocity * delta_t);

        p(planetA,:) = newPlanetAPosition;
        v(planetA,:) = newPlanetAVelocity;
        disp(newPlanetAPosition);
        % Draw positions of planet
        addpoints(pens(planetA), p(planetA,1), p(planetA,2), p(planetA,3));
        drawnow(); 
    end
    
    
    
end
      
end
    

   
    