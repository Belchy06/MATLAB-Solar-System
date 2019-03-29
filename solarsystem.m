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
% Graph Variables
pens = [animatedline('LineWidth', 2),animatedline('LineWidth', 2),animatedline('LineWidth', 2),animatedline('LineWidth', 2),animatedline('LineWidth', 2)];
set(gca,'XLim', [-1836100000000 1836100000000],'YLim', [-1836100000000 1836100000000],'ZLim', [-1836100000000 1836100000000]);  
 

% Time constants
delta_t = 60 * 60 * 24;

% Universal Gravitational Constant
G = 6.673 * 10^-11;
    
% Graph Variables
 
% Write your code here
% Loop through code for each second
for t=1:stop_time
   for planetA=1:length(p)
        planetAPosition = [p(planetA,1), p(planetA,2)];
        planetAVelocity = [v(planetA,1), v(planetA,2)];
        planetAMass = mass(planetA);
        
       for planetB=1:length(p) 
            planetBPosition = [p(planetB,1), p(planetB,2)];
            planetBMass = mass(planetB);

            if planetA ~= planetB
                differenceVector = planetBPosition - planetAPosition;
                differenceMagnitude = norm(differenceVector);
                
                forceVector = G * ((planetAMass * planetBMass) / differenceMagnitude^3) .* differenceVector;
                
                % Gravity Vector from planet B -> A
                % Positive value as specified in assginment sheet
                gravityVectorArray.x(planetB,planetA) = forceVector(1);
                gravityVectorArray.y(planetB,planetA) = forceVector(2);

                % Gravity Vector from planet A -> B
                % Negative value can be achieved by multiplying gravity vector component by -1
                % as the direction is just reversed. No need to recalculate components
                gravityVectorArray.x(planetA,planetB) = forceVector(1) * -1;
                gravityVectorArray.y(planetA,planetB) = forceVector(2) * -1;
                
                
                newAcceleration.x = (1 / planetAMass) * gravityVectorArray(:,planetA).x;
                newAcceleration.y = (1 / planetAMass) * gravityVectorArray(:,planetA).y;
                
                newVelocity.x = planetAVelocity(1) + (newAcceleration.x * delta_t);
                newVelocity.y = planetAVelocity(2) + (newAcceleration.y * delta_t);
                
                newPosition.x = (planetAVelocity(1) * delta_t) + (0.5 * newAcceleration.x * (delta_t^2));
                newPosition.y = (planetAVelocity(2) * delta_t) + (0.5 * newAcceleration.y * (delta_t^2));
                
                p(planetA) = newPosition;
                v(planetA) = newVelocity;
                
                
                addpoints(pens(planetA), p(planetA,1), p(planetA,2), p(planetA,3));
                drawnow(); 
            end      
       end
   end
end
  

end
    

   
    