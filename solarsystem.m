function [p, v, e] = solarsystem(p, v, e, mass, stop_time, hide_animation)
  
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
    delta_t = 60*60;
    
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
        for planetI=1:length(p)
            planetIMass = mass(planetI);
            planetIPosition = [p(planetI,1), p(planetI,2)];
            planetIVelocity = [v(planetI,1), v(planetI,2)];
    
            for planetJ=(planetI+1):length(p)
                planetJMass = mass(planetJ)
                planetJPosition = [p(planetJ,1), p(planetJ,2)];
                distanceVector = planetJPosition - planetIPosition;
                distanceVectorMagnitude = sqrt(distanceVector(1)^2 + distanceVector(2)^2);
    
                % Calculate Potential Energy between two objects
                potentialEnergy = -G * ((planetIMass * planetJMass) / distanceVectorMagnitude);
                initialEnergyArray(planetJ, planetI) = potentialEnergy;
                initialEnergyArray(planetI, planetJ) = potentialEnergy;
            end
    
            kineticEnergy(planetI) = planetIKinetic;
            potentialEnergy(planetI) = sum(initialEnergyArray(:,planetI));
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
        for planetI=1:length(p)
            planetIMass = mass(planetI);
            planetIPosition = [p(planetI,1), p(planetI,2), p(planetI,3)];
            planetIVelocity = [v(planetI,1), v(planetI,2), v(planetI,3)];
            planetIKinetic = 0.5 * planetIMass * sqrt(planetIVelocity(1)^2 + planetIVelocity(2)^2 + planetIVelocity(3)^2)^2;
    
            for planetJ=(planetI+1):length(p)
                planetJMass = mass(planetJ);
                planetJPosition = [p(planetJ,1), p(planetJ,2), p(planetJ,3)];
                distanceVector = planetJPosition - planetIPosition;
                distanceVectorMagnitude = sqrt(distanceVector(1)^2 + distanceVector(2)^2 + distanceVector(3)^2);
    
                % Calculate Potential Energy between two objects
                potentialEnergy = -G * ((planetIMass * planetJMass) / distanceVectorMagnitude);
                initialEnergyArray(planetJ, planetI) = potentialEnergy;
                initialEnergyArray(planetI, planetJ) = potentialEnergy;
            end
    
            kineticEnergy(planetI) = planetIKinetic;
            potentialEnergy(planetI) = sum(initialEnergyArray(:,planetI));
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
        for t=1:delta_t:stop_time   
            for planetI=1:length(p)
                % Retrieve data from supplied data
                planetIPosition = [p(planetI,1), p(planetI,2)];
                planetIVelocity = [v(planetI,1), v(planetI,2)];
                planetIMass = mass(planetI);
                
                for planetJ=(planetI+1):length(p)      
                    % Retrieve data from supplied data
                    planetJPosition = [p(planetJ,1), p(planetJ,2)];
                    planetJMass = mass(planetJ);
    
                    % Calculate necessary values as defined in assignment sheet
                    distanceVector = planetJPosition - planetIPosition;
                    distanceVectorMagnitude = sqrt(distanceVector(1)^2 + distanceVector(2)^2);
    
                    gravityVector = ((G * planetIMass * planetJMass) / (distanceVectorMagnitude^3)) * distanceVector;
                
                    % Gravity Vector from planet B -> A
                    gravityVectorArray.x(planetJ,planetI) = gravityVector(1);
                    gravityVectorArray.y(planetJ,planetI) = gravityVector(2);
    
                    % Gravity Vector from planet A -> B is just the inverse of
                    % B -> A
                    gravityVectorArray.x(planetI,planetJ) = gravityVector(1) * -1;
                    gravityVectorArray.y(planetI,planetJ) = gravityVector(2) * -1;
                end
    
                % Total Gravity in xyz components that affects planetI
                totalGravity.x = sum(gravityVectorArray.x(:,planetI));
                totalGravity.y = sum(gravityVectorArray.y(:,planetI));
            
                % Move total gravity into an array of acceleration values
                gravitationalAcceleration = [totalGravity.x/planetIMass totalGravity.y/planetIMass];
                newplanetIVelocity = planetIVelocity + (delta_t * gravitationalAcceleration);
                newplanetIPosition = planetIPosition + (newplanetIVelocity * delta_t) + (gravitationalAcceleration .* (0.5 * delta_t^2));
    
                p(planetI,:) = newplanetIPosition;
                v(planetI,:) = newplanetIVelocity;
    
                % Draw positions of planet
                plot2DPlanet(pens(planetI), p(planetI, 1), p(planetI, 2)); 
            end   
        end
        calculateTotal2DEnergy();
    end
    
    function animate3DSystem()
        % Loop through code for each second
        for t=1:delta_t:stop_time   
            for planetI=1:length(p)    
                % Retrieve data from supplied data
                planetIPosition = [p(planetI,1), p(planetI,2), p(planetI,3)];
                planetIVelocity = [v(planetI,1), v(planetI,2), v(planetI,3)];
                planetIMass = mass(planetI);
                
                for planetJ=(planetI+1):length(p)
                    % Retrieve data from supplied data
                    planetJPosition = [p(planetJ,1), p(planetJ,2), p(planetJ,3)];
                    planetJMass = mass(planetJ);
    
                    % Calculate necessary values as defined in assignment sheet
                    distanceVector = planetJPosition - planetIPosition;
                    distanceVectorMagnitude = sqrt(distanceVector(1)^2 + distanceVector(2)^2 + distanceVector(3)^2);
    
                    gravityVector = ((G * planetIMass * planetJMass) / (distanceVectorMagnitude^3)) * distanceVector;
                
                    % Gravity Vector from planet B -> A
                    % Positive value as specified in assginment sheet
                    gravityVectorArray.x(planetJ,planetI) = gravityVector(1);
                    gravityVectorArray.y(planetJ,planetI) = gravityVector(2);
                    gravityVectorArray.z(planetJ,planetI) = gravityVector(3);
    
                    % Gravity Vector from planet A -> B
                    % Negative value can be achieved by multiplying gravity vector component by -1
                    % as the direction is just reversed. No need to recalculate components
                    gravityVectorArray.x(planetI,planetJ) = gravityVector(1) * -1;
                    gravityVectorArray.y(planetI,planetJ) = gravityVector(2) * -1;
                    gravityVectorArray.z(planetI,planetJ) = gravityVector(3) * -1;                 
                end
                % Total Gravity in xyz components that affects planetI
                totalGravity.x = sum(gravityVectorArray.x(:,planetI));
                totalGravity.y = sum(gravityVectorArray.y(:,planetI));
                totalGravity.z = sum(gravityVectorArray.z(:,planetI));
                
                % Move total gravity into an array of acceleration values
                gravitationalAcceleration = [totalGravity.x/planetIMass totalGravity.y/planetIMass totalGravity.z/planetIMass];
                newplanetIVelocity = planetIVelocity + (delta_t * gravitationalAcceleration);
                newplanetIPosition = planetIPosition + (newplanetIVelocity * delta_t) + (gravitationalAcceleration .* (0.5 * delta_t^2));
    
                p(planetI,:) = newplanetIPosition;
                v(planetI,:) = newplanetIVelocity;
               
                % Draw positions of planet
                plot3DPlanet(pens(planetI), p(planetI, 1), p(planetI, 2), p(planetI, 3));     
            end   
        end
        calculateTotal3DEnergy();
    end
    
    
    % Simulate Solar System
    function simulate2DSystem()
        % Loop through code for each second
        for t=1:delta_t:stop_time   
            for planetI=1:length(p)
                % Retrieve data from supplied data
                planetIPosition = [p(planetI,1), p(planetI,2)];
                planetIVelocity = [v(planetI,1), v(planetI,2)];
                planetIMass = mass(planetI);
    
                for planetJ=(planetI+1):length(p)      
                    % Retrieve data from supplied data
                    planetJPosition = [p(planetJ,1), p(planetJ,2)];
                    planetJMass = mass(planetJ);
    
                    % Calculate necessary values as defined in assignment sheet
                    distanceVector = planetJPosition - planetIPosition;
                    distanceVectorMagnitude = sqrt(distanceVector(1)^2 + distanceVector(2)^2);
    
                    gravityVector = ((G * planetIMass * planetJMass) / (distanceVectorMagnitude^3)) * distanceVector;
                
                    % Gravity Vector from planet B -> A
                    gravityVectorArray.x(planetJ,planetI) = gravityVector(1);
                    gravityVectorArray.y(planetJ,planetI) = gravityVector(2);
    
                    % Gravity Vector from planet A -> B is just the inverse of
                    % B -> A
                    gravityVectorArray.x(planetI,planetJ) = gravityVector(1) * -1;
                    gravityVectorArray.y(planetI,planetJ) = gravityVector(2) * -1;
                end
    
                % Total Gravity in xyz components that affects planetI
                totalGravity.x = sum(gravityVectorArray.x(:,planetI));
                totalGravity.y = sum(gravityVectorArray.y(:,planetI));
            
                % Move total gravity into an array of acceleration values
                gravitationalAcceleration = [totalGravity.x/planetIMass totalGravity.y/planetIMass];
                newplanetIVelocity = planetIVelocity + (delta_t * gravitationalAcceleration);
                newplanetIPosition = planetIPosition + (newplanetIVelocity * delta_t) + (gravitationalAcceleration .* (0.5 * delta_t^2));
    
                p(planetI,:) = newplanetIPosition;
                v(planetI,:) = newplanetIVelocity;
            end   
        end
        calculateTotal2DEnergy();
    end
    
    function simulate3DSystem()
        % Loop through code for each second
        for t=1:delta_t:stop_time   
            for planetI=1:length(p)    
                % Retrieve data from supplied data
                planetIPosition = [p(planetI,1), p(planetI,2), p(planetI,3)];
                planetIVelocity = [v(planetI,1), v(planetI,2), v(planetI,3)];
                planetIMass = mass(planetI);
    
                for planetJ=(planetI+1):length(p)
                    % Retrieve data from supplied data
                    planetJPosition = [p(planetJ,1), p(planetJ,2), p(planetJ,3)];
                    planetJMass = mass(planetJ);
    
                    % Calculate necessary values as defined in assignment sheet
                    distanceVector = planetJPosition - planetIPosition;
                    distanceVectorMagnitude = sqrt(distanceVector(1)^2 + distanceVector(2)^2 + distanceVector(3)^2);
    
                    gravityVector = ((G * planetIMass * planetJMass) / (distanceVectorMagnitude^3)) * distanceVector;
                
                    % Gravity Vector from planet B -> A
                    % Positive value as specified in assginment sheet
                    gravityVectorArray.x(planetJ,planetI) = gravityVector(1);
                    gravityVectorArray.y(planetJ,planetI) = gravityVector(2);
                    gravityVectorArray.z(planetJ,planetI) = gravityVector(3);
    
                    % Gravity Vector from planet A -> B
                    % Negative value can be achieved by multiplying gravity vector component by -1
                    % as the direction is just reversed. No need to recalculate components
                    gravityVectorArray.x(planetI,planetJ) = gravityVector(1) * -1;
                    gravityVectorArray.y(planetI,planetJ) = gravityVector(2) * -1;
                    gravityVectorArray.z(planetI,planetJ) = gravityVector(3) * -1;                   
                end
    
                % Total Gravity in xyz components that affects planetI
                totalGravity.x = sum(gravityVectorArray.x(:,planetI));
                totalGravity.y = sum(gravityVectorArray.y(:,planetI));
                totalGravity.z = sum(gravityVectorArray.z(:,planetI));
            
                % Move total gravity into an array of acceleration values
                gravitationalAcceleration = [totalGravity.x/planetIMass totalGravity.y/planetIMass totalGravity.z/planetIMass];
                newplanetIVelocity = planetIVelocity + (delta_t * gravitationalAcceleration);
                newplanetIPosition = planetIPosition + (newplanetIVelocity * delta_t) + (gravitationalAcceleration .* (0.5 * delta_t^2));
    
                p(planetI,:) = newplanetIPosition;
                v(planetI,:) = newplanetIVelocity;
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
        
    
       
        