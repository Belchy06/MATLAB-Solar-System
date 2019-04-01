function [e] = solarsystem(p, v, e, mass, stop_time, hide_animation)
  
    if nargin < 6
        hide_animation = false;
    end
    
    % Graph Constants
    graphBounds = 3 * 10^11;
    set(gcf, 'name', 'Solar System Model');
    xlabel('X (m)'); 
    ylabel('Y (m)'); 
    zlabel('Z (m)');
    set(gca,'XLim', [-graphBounds graphBounds],'YLim', [-graphBounds graphBounds],'ZLim', [-graphBounds graphBounds]);  
    pens = [animatedline('LineWidth', 2),
            animatedline('LineWidth', 4, 'Color','c'),
            animatedline('LineWidth', 2, 'Color',[1 0.5 0]),
            animatedline('LineWidth', 3, 'Color','y'),
            animatedline('LineWidth', 4, 'Color','r'),
            animatedline('LineWidth', 2),
            animatedline('LineWidth', 2),
            animatedline('LineWidth', 2),
            animatedline('LineWidth', 2),
            animatedline('LineWidth', 2),
            animatedline('LineWidth', 4)];
    
    % Time Constant
    delta_t = 1;
    
    % Universal Gravitational Constant
    G = 6.673 * 10^-11;
    
    % Initialize Variables
    initialEnergy = 0;
    finalEnergy = 0;
    kineticEnergy = 0;
    potentialEnergy = 0;
    firstRun = true;
    
    
    % Start of logic
    ArrayDimensions = size(p,2);
    switch ArrayDimensions
        case 2
            calculateTotal2DEnergy();
            if hide_animation == false
                animate2DSystem();
            else
                simulate2DSystem();
            end        
            
        case 3
            calculateTotal3DEnergy();
            if hide_animation == false
                animate3DSystem();
            else
                simulate3DSystem();
            end
            
            
        otherwise
            disp('Expected input dimensions of "p" to be %ix2 or %ix3, received %ix%i instead.', size(p,1), size(p,1), size(p,1), size(p,2));
    end
    
    function calculateTotal2DEnergy()
        for planetA=1:length(p)
            planetAMass = mass(planetA);
            planetAPosition = [p(planetA,1), p(planetA,2)];
            planetAVelocity = [v(planetA,1), v(planetA,2)];
    
            for planetB=(planetA+1):length(p)
                planetBMass = mass(planetB)
                planetBPosition = [p(planetB,1), p(planetB,2)];
                distanceVector = planetBPosition - planetAPosition;
                distanceVectorMagnitude = sqrt(distanceVector(1)^2 + distanceVector(2)^2);
    
                % Calculate Potential Energy between two objects
                potentialEnergy = -G * ((planetAMass * planetBMass) / distanceVectorMagnitude);
                initialEnergyArray(planetB, planetA) = potentialEnergy;
                initialEnergyArray(planetA, planetB) = potentialEnergy;
            end
    
            kineticEnergy(planetA) = planetAKinetic;
            potentialEnergy(planetA) = sum(initialEnergyArray(:,planetA));
        end
    
        if firstRun
            initialEnergy = sum(kineticEnergy) + sum(potentialEnergy);
            kineticEnergy = 0;
            potentialEnergy = 0;
            firstRun = false;
        else
            finalEnergy = sum(kineticEnergy) + sum(potentialEnergy);
            calculateEnergyError();
        end
    end
    
    function calculateTotal3DEnergy()
        for planetA=1:length(p)
            planetAMass = mass(planetA);
            planetAPosition = [p(planetA,1), p(planetA,2), p(planetA,3)];
            planetAVelocity = [v(planetA,1), v(planetA,2), v(planetA,3)];
            planetAKinetic = 0.5 * planetAMass * sqrt(planetAVelocity(1)^2 + planetAVelocity(2)^2 + planetAVelocity(3)^2)^2;
    
            for planetB=(planetA+1):length(p)
                planetBMass = mass(planetB);
                planetBPosition = [p(planetB,1), p(planetB,2), p(planetB,3)];
                distanceVector = planetBPosition - planetAPosition;
                distanceVectorMagnitude = sqrt(distanceVector(1)^2 + distanceVector(2)^2 + distanceVector(3)^2);
    
                % Calculate Potential Energy between two objects
                potentialEnergy = -G * ((planetAMass * planetBMass) / distanceVectorMagnitude);
                initialEnergyArray(planetB, planetA) = potentialEnergy;
                initialEnergyArray(planetA, planetB) = potentialEnergy;
            end
    
            kineticEnergy(planetA) = planetAKinetic;
            potentialEnergy(planetA) = sum(initialEnergyArray(:,planetA));
        end
    
        if firstRun == true
            initialEnergy = sum(kineticEnergy) + sum(potentialEnergy);
            kineticEnergy = 0;
            potentialEnergy = 0;
            firstRun = false;
        else
            finalEnergy = sum(kineticEnergy) + sum(potentialEnergy);
            calculateEnergyError();
        end
    end
    
    % Animate Solar System
    function animate2DSystem()
        % Loop through code for each second
        for t=1:stop_time   
            for planetA=1:length(p)
                % Retrieve data from supplied data
                planetAPosition = [p(planetA,1), p(planetA,2)];
                planetAVelocity = [v(planetA,1), v(planetA,2)];
                planetAMass = mass(planetA);
                
                for planetB=(planetA+1):length(p)      
                    % Retrieve data from supplied data
                    planetBPosition = [p(planetB,1), p(planetB,2)];
                    planetBMass = mass(planetB);
    
                    % Calculate necessary values as defined in assignment sheet
                    distanceVector = planetBPosition - planetAPosition;
                    distanceVectorMagnitude = sqrt(distanceVector(1)^2 + distanceVector(2)^2);
    
                    gravityVector = ((G * planetAMass * planetBMass) / (distanceVectorMagnitude^3)) * distanceVector;
                
                    % Gravity Vector from planet B -> A
                    gravityVectorArray.x(planetB,planetA) = gravityVector(1);
                    gravityVectorArray.y(planetB,planetA) = gravityVector(2);
    
                    % Gravity Vector from planet A -> B is just the inverse of
                    % B -> A
                    gravityVectorArray.x(planetA,planetB) = gravityVector(1) * -1;
                    gravityVectorArray.y(planetA,planetB) = gravityVector(2) * -1;
                end
    
                % Total Gravity in xyz components that affects planetA
                totalGravity.x = sum(gravityVectorArray.x(:,planetA));
                totalGravity.y = sum(gravityVectorArray.y(:,planetA));
            
                % Move total gravity into an array of acceleration values
                gravitationalAcceleration = [totalGravity.x/planetAMass totalGravity.y/planetAMass];
                newPlanetAVelocity = planetAVelocity + (delta_t * gravitationalAcceleration);
                newPlanetAPosition = planetAPosition + (newPlanetAVelocity * delta_t) + (gravitationalAcceleration .* (0.5 * delta_t^2));
    
                p(planetA,:) = newPlanetAPosition;
                v(planetA,:) = newPlanetAVelocity;
    
                % Draw positions of planet
                plot2DPlanet(pens(planetA), p(planetA, 1), p(planetA, 2)); 
            end   
        end
        calculateTotal2DEnergy();
    end
    
    function animate3DSystem()
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
    
                    gravityVector = ((G * planetAMass * planetBMass) / (distanceVectorMagnitude^3)) * distanceVector;
                
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
                newPlanetAVelocity = planetAVelocity + (delta_t * gravitationalAcceleration);
                newPlanetAPosition = planetAPosition + (newPlanetAVelocity * delta_t) + (gravitationalAcceleration .* (0.5 * delta_t^2));
    
                p(planetA,:) = newPlanetAPosition;
                v(planetA,:) = newPlanetAVelocity;
               
                % Draw positions of planet
                plot3DPlanet(pens(planetA), p(planetA, 1), p(planetA, 2), p(planetA, 3));     
            end   
        end
        calculateTotal3DEnergy();
    end
    
    
    % Simulate Solar System
    function simulate2DSystem()
        % Loop through code for each second
        for t=1:stop_time   
            for planetA=1:length(p)
                % Retrieve data from supplied data
                planetAPosition = [p(planetA,1), p(planetA,2)];
                planetAVelocity = [v(planetA,1), v(planetA,2)];
                planetAMass = mass(planetA);
    
                for planetB=(planetA+1):length(p)      
                    % Retrieve data from supplied data
                    planetBPosition = [p(planetB,1), p(planetB,2)];
                    planetBMass = mass(planetB);
    
                    % Calculate necessary values as defined in assignment sheet
                    distanceVector = planetBPosition - planetAPosition;
                    distanceVectorMagnitude = sqrt(distanceVector(1)^2 + distanceVector(2)^2);
    
                    gravityVector = ((G * planetAMass * planetBMass) / (distanceVectorMagnitude^3)) * distanceVector;
                
                    % Gravity Vector from planet B -> A
                    gravityVectorArray.x(planetB,planetA) = gravityVector(1);
                    gravityVectorArray.y(planetB,planetA) = gravityVector(2);
    
                    % Gravity Vector from planet A -> B is just the inverse of
                    % B -> A
                    gravityVectorArray.x(planetA,planetB) = gravityVector(1) * -1;
                    gravityVectorArray.y(planetA,planetB) = gravityVector(2) * -1;
                end
    
                % Total Gravity in xyz components that affects planetA
                totalGravity.x = sum(gravityVectorArray.x(:,planetA));
                totalGravity.y = sum(gravityVectorArray.y(:,planetA));
            
                % Move total gravity into an array of acceleration values
                gravitationalAcceleration = [totalGravity.x/planetAMass totalGravity.y/planetAMass];
                newPlanetAVelocity = planetAVelocity + (delta_t * gravitationalAcceleration);
                newPlanetAPosition = planetAPosition + (newPlanetAVelocity * delta_t) + (gravitationalAcceleration .* (0.5 * delta_t^2));
    
                p(planetA,:) = newPlanetAPosition;
                v(planetA,:) = newPlanetAVelocity;
            end   
        end
        calculateTotal2DEnergy();
    end
    
    function simulate3DSystem()
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
    
                    gravityVector = ((G * planetAMass * planetBMass) / (distanceVectorMagnitude^3)) * distanceVector;
                
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
                newPlanetAVelocity = planetAVelocity + (delta_t * gravitationalAcceleration);
                newPlanetAPosition = planetAPosition + (newPlanetAVelocity * delta_t) + (gravitationalAcceleration .* (0.5 * delta_t^2));
    
                p(planetA,:) = newPlanetAPosition;
                v(planetA,:) = newPlanetAVelocity;
            end   
        end
        calculateTotal3DEnergy();
    end
    
    % Plot Planets
    function plot2DPlanet(pen, x, y)
        addpoints(pen, x, y);     
        drawnow(); 
    end
    
    function plot3DPlanet(pen, x, y, z)
        addpoints(pen, x, y, z);     
        drawnow(); 
    end
    
    % Calculate percentage error between inital energy in system and final
    % energy in system
    function calculateEnergyError()
        % Calculate Final Energy Values
        e = abs(finalEnergy - initialEnergy) / initialEnergy;
    end
    
    end
        
    
       
        