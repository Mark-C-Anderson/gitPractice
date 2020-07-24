# Purpose: Model a geometry using a Vortex Lattice Method
# Inputs:
#   panels: the panels used to create the geometry of the wing. Most easily generated properly
#           by using the generatePanels() function
#   freestream: the three-dimensional freestream flow at each panel
#   density: the ambient fluid density
#
# Outputs:
#   CL: total lift coefficient
#   CDi_near: induced drag coefficient as calculated using Kutta-Joukowski in the near field
#   cl: sectional lift Coefficients
#   cd: sectional induced drag Coefficients
#   CLSpanLocations: the y-coordinate for the center of each panel
#   GammaValues: The circulation values for each panel

# To use this code, you need to be able to run julia files.

# This is the main module
module VLMMCA

    using Pkg
    Pkg.add("LinearAlgebra")
    Pkg.add("PyPlot")
    using LinearAlgebra
    using PyPlot

    function VLM(panels,
                 freestream,
                 density = 1.225)

        # Preparing the AIC matrix
        AIC,unitNormals = createAIC(panels,"Ring Lattice")

        # Defining the b column vector
        b = zeros(length(unitNormals[:,1]))
        for i = 1:length(b)
            u = freestream[i,:]
            n = [unitNormals[i,1],unitNormals[i,2],unitNormals[i,3]]
            b[i] = -dot(u,n)
        end
        
        # Solving for the circulations
        GammaValues = AIC\b

        # Calculating the lift
        CL, cl, CLSpanLocations = calculateLift(density,freestream,panels,GammaValues);

        # Calculating the induced drag in the near field
        CDi_near, cd_near, Cd_nearSpanLocations, inducedVelocity = calculateInducedDrag(density,freestream,panels,GammaValues,cl);

        # Calculating the induced drag in the far field
        CDi_far, cd_far, Cd_farSpanLocations = calculateDrag(density,freestream,panels,GammaValues);

        return CL, CDi_near, cl, cd_near, CLSpanLocations, GammaValues

    end

    function VelocityFilament(r1,r2,gamma,type)

        if type == "Bound"
    
            V = cross(r1,r2)./(norm(r1)*norm(r2) + dot(r1,r2))
            V = V .* (1/norm(r1)+1/norm(r2))
            V = V.*gamma./(4 * pi)
    
        elseif type == "Left"
    
            V = cross(r1,[1,0,0])./(norm(r1) - dot(r1,[1,0,0]))
            V = V .* (1/norm(r1))
            V = V.*gamma./(4 * pi)
    
        elseif type == "Right"
            #println("r2 = ", r2)
            V = cross(r2,[1,0,0])./(norm(r2) - dot(r2,[1,0,0]))
            V = V .* (1/norm(r2))
            V = -V.*gamma./(4 * pi)
    
        end
    
        return V
    
    end

    function definePoints(panels)

        coordinates = zeros(length(panels[:,1]),9)
        
        for i = 1:length(panels[:,1])
            currentPanel = panels[i,:]
        
            # To calculate bound vortex x-coordinates:
            # Average of two x-coordinates - 1/4* Difference of the same two x-coordinates
            X1n = (currentPanel[1] + currentPanel[10])/2 - (1/4)*abs(currentPanel[1] - currentPanel[10])
            X2n = (currentPanel[4] + currentPanel[7])/2 - (1/4)*abs(currentPanel[4] - currentPanel[7])
        
            # To calculate bound vortex y-coordinates:
            Y1n = (currentPanel[2] + currentPanel[11])/2 - (1/4)*abs(currentPanel[2] - currentPanel[11])
            Y2n = (currentPanel[5] + currentPanel[8])/2 - (1/4)*abs(currentPanel[5] - currentPanel[8])
        
            # To calculate the bound vortex z-coordinates:
            Z1n = (currentPanel[3] + currentPanel[12])/2 - (1/4)*abs(currentPanel[3] - currentPanel[12])
            Z2n = (currentPanel[6] + currentPanel[9])/2 - (1/4)*abs(currentPanel[6] - currentPanel[9]);
            
            # To calculate control point coordinates:
            # Use the same method as in finding the bound vortex x-coordinates, but this time
            # do both sides of the panel and average their values, giving a point along the axis
            # of the panel.
            # Do the same for the y-coordinate
            # This part is important for the airplane geometry
            Xm = ((currentPanel[1] + currentPanel[10])/2 + (1/4)*abs(currentPanel[1] - currentPanel[10]) + 
                (currentPanel[4] + currentPanel[7])/2 + (1/4)*abs(currentPanel[4] - currentPanel[7]))/2
            Ym = ((currentPanel[2] + currentPanel[11])/2 + (1/4)*abs(currentPanel[2] - currentPanel[11]) + 
                (currentPanel[5] + currentPanel[8])/2 + (1/4)*abs(currentPanel[5] - currentPanel[8]))/2
            Zm = ((currentPanel[3] + currentPanel[12])/2 + (1/4)*abs(currentPanel[3] - currentPanel[12]) + 
                (currentPanel[6] + currentPanel[9])/2 + (1/4)*abs(currentPanel[6] - currentPanel[9]))/2
            
            coordinates[i,:] = [Xm,Ym,Zm,X1n,Y1n,Z1n,X2n,Y2n,Z2n]

        end
        
        coordinates
        
    end

    function Velocity(r1,r2,gammaValue,vortexType)

        if vortexType == "Horseshoe"
            velocity = cross(r1,r2)./(norm(r1)*norm(r2) + dot(r1,r2))
            velocity = velocity .* (1/norm(r1)+1/norm(r2))
            velocity = velocity + cross(r1,[1,0,0])./(norm(r1)-r1[1]).*1/norm(r1)
            velocity = velocity - cross(r2,[1,0,0])./(norm(r2)-r2[1]).*1/norm(r2)
    
            velocity = velocity * gammaValue / (4*pi)
    
        end
    
        if vortexType == "Bound"
            velocity = cross(r1,r2)./(norm(r1)*norm(r2) + dot(r1,r2))
            velocity = velocity .* (1/norm(r1)+1/norm(r2))
            velocity = velocity.*gammaValue./(4 * pi)
        end
    
        if vortexType == "Semi-Infinite Left"
            velocity = cross(r1,[1,0,0])./(norm(r1) - dot(r1,[1,0,0]))
            velocity = velocity .* (1/norm(r1))
            velocity = velocity.*gammaValue./(4 * pi)
        end
    
        if vortexType == "Semi-Infinite Right"
            velocity = cross(r2,[1,0,0])./(norm(r2) - dot(r2,[1,0,0]))
            velocity = velocity .* (1/norm(r2))
            velocity = -velocity.*gammaValue./(4 * pi)
        end
    
        return velocity
    
    end

    function generatePanels(firstCoordinate, secondCoordinate, thirdCoordinate, fourthCoordinate, numPanelsSpan, numPanelsChord = 1)
        # Coordinates must be defined starting with the point on the leading edge nearest the center of the wing, and then
        # proceed clockwise as viewed looking down on the wing.
        # The variable "numPanels" is the number of panels on just half of the total wing, or in other words is the number
        # of panels across the coordiantes you define.
        
        # Initializing the finalPanels array that will be the returned value of this function
        finalPanels = []

        # Breaking the coordinates into chordwise values
        rootXSpacing = (fourthCoordinate[1] - firstCoordinate[1])/numPanelsChord
        tipXSpacing  = (thirdCoordinate[1] - secondCoordinate[1])/numPanelsChord

        rootYSpacing = (fourthCoordinate[2] - firstCoordinate[2])/numPanelsChord
        tipYSpacing  = (thirdCoordinate[2] - secondCoordinate[2])/numPanelsChord

        rootZSpacing = (fourthCoordinate[3] - firstCoordinate[3])/numPanelsChord
        tipZSpacing  = (thirdCoordinate[3] - secondCoordinate[3])/numPanelsChord

        # Creating arrays to store the different coordinate values in
        adjustedFirstCoordinate = zeros(numPanelsChord,3)
        adjustedSecondCoordinate = zeros(numPanelsChord,3)
        adjustedThirdCoordinate = zeros(numPanelsChord,3)
        adjustedFourthCoordinate = zeros(numPanelsChord,3)

        # Assigning the new coordinate values to the arrays to account for chordwise panels
        for i = 1:numPanelsChord

            # We want the first and second coordinates to march toward the trailing edge of the wing
            adjustedFirstCoordinate[i,1] = firstCoordinate[1] + (i-1)*rootXSpacing
            adjustedFirstCoordinate[i,2] = firstCoordinate[2] + (i-1)*rootYSpacing
            adjustedFirstCoordinate[i,3] = firstCoordinate[3] + (i-1)*rootZSpacing

            adjustedSecondCoordinate[i,1] = secondCoordinate[1] + (i-1)*tipXSpacing
            adjustedSecondCoordinate[i,2] = secondCoordinate[2] + (i-1)*tipYSpacing
            adjustedSecondCoordinate[i,3] = secondCoordinate[3] + (i-1)*tipZSpacing

            # We want the third and fourth coordinates to march toward the trailing edge too, but to start closer to the leading edge
            adjustedThirdCoordinate[i,1] = adjustedSecondCoordinate[i,1] + tipXSpacing
            adjustedThirdCoordinate[i,2] = adjustedSecondCoordinate[i,2] + tipYSpacing
            adjustedThirdCoordinate[i,3] = adjustedSecondCoordinate[i,3] + tipZSpacing

            adjustedFourthCoordinate[i,1] = adjustedFirstCoordinate[i,1] + rootXSpacing
            adjustedFourthCoordinate[i,2] = adjustedFirstCoordinate[i,2] + rootYSpacing
            adjustedFourthCoordinate[i,3] = adjustedFirstCoordinate[i,3] + rootZSpacing
            

        end

        # println(adjustedFirstCoordinate)
        # println(adjustedSecondCoordinate)
        # println(adjustedThirdCoordinate)
        # println(adjustedFourthCoordinate)

        # Reassigning the coordinates to work with the old code
        firstCoordinate = adjustedFirstCoordinate
        secondCoordinate = adjustedSecondCoordinate
        thirdCoordinate = adjustedThirdCoordinate
        fourthCoordinate = adjustedFourthCoordinate

        for j = 1:numPanelsChord
            # initializing the final panels array
            panels = zeros(2*numPanelsSpan,13)
            orderedPanels = zeros(2*numPanelsSpan,13)
            
            # Getting the various spacings sorted out
            leadingEdgeXSpacing = (secondCoordinate[j,1] - firstCoordinate[j,1])/numPanelsSpan;
            trailingEdgeXSpacing = (thirdCoordinate[j,1] - fourthCoordinate[j,1])/numPanelsSpan;
            
            leadingEdgeYSpacing = (secondCoordinate[j,2] - firstCoordinate[j,2])/numPanelsSpan;
            trailingEdgeYSpacing = (thirdCoordinate[j,2] - fourthCoordinate[j,2])/numPanelsSpan;
            
            leadingEdgeZSpacing = (secondCoordinate[j,3] - firstCoordinate[j,3])/numPanelsSpan;
            trailingEdgeZSpacing = (thirdCoordinate[j,3] - fourthCoordinate[j,3])/numPanelsSpan;
            
            
            # Defining the x-coordinates on the right side of the wing that the coordinates were given for
            for i = 1:numPanelsSpan
                    
                panels[i,1] = firstCoordinate[j,1] + (i-1)*leadingEdgeXSpacing; # Top left corner
                panels[i,4] = panels[i,1] + leadingEdgeXSpacing; # Top right corner
                panels[i,10] = fourthCoordinate[j,1] + (i-1)*trailingEdgeXSpacing; # Bottom left corner
                panels[i,7] = panels[i,10] + trailingEdgeXSpacing; # Bottom left corner

                # Setting a flag to tell us whether this is a trailing edge panel
                if j == numPanelsChord
                    panels[i,13] = 1
                end
                
            end
            
            # Defining the y-coordinates on the right side of the wing that the coordinates were given for
            for i = 1:numPanelsSpan
                    
                panels[i,2] = firstCoordinate[j,2] + (i-1)*leadingEdgeYSpacing; # Top left corner
                panels[i,5] = panels[i,2] + leadingEdgeYSpacing; # Top right corner
                panels[i,11] = fourthCoordinate[j,2] + (i-1)*trailingEdgeYSpacing; # Bottom left corner
                panels[i,8] = panels[i,11] + trailingEdgeYSpacing; # Bottom left corner

                # Setting a flag to tell us whether this is a trailing edge panel
                if j == numPanelsChord
                    panels[i,13] = 1
                end
                
            end
            
            # Defining the z-coordinates on the right side of the wing that the coordinates were given for
            for i = 1:numPanelsSpan
                    
                panels[i,3] = firstCoordinate[j,3] + (i-1)*leadingEdgeZSpacing; # Top left corner
                panels[i,6] = panels[i,3] + leadingEdgeZSpacing; # Top right corner
                panels[i,12] = fourthCoordinate[j,3] + (i-1)*trailingEdgeZSpacing; # Bottom left corner
                panels[i,9] = panels[i,12] + trailingEdgeZSpacing; # Bottom left corner

                # Setting a flag to tell us whether this is a trailing edge panel
                if j == numPanelsChord
                    panels[i,13] = 1
                end
                
            end
            
            # Defining the symmetric side of wing
            for i = 1:numPanelsSpan
                
                # The wing is assumed to be symmetric about the x-axis, such that the y-values must be negative.
                # In order to maintain consistent orientation of the unit normal to each panel, the order of the
                # coordinates for each panel must be flipped such that the first and second, along with the third
                # and fourth swap places.
                
                # Top Left (swapping first and second, and the y-value becomes negative)
                panels[numPanelsSpan + i,1] = panels[i,4]
                panels[numPanelsSpan + i,2] = -panels[i,5]
                panels[numPanelsSpan + i,3] = panels[i,6]
                
                # Top Right (swapping first and second, and the y-value becomes negative)
                panels[numPanelsSpan + i,4] = panels[i,1]
                panels[numPanelsSpan + i,5] = -panels[i,2]
                panels[numPanelsSpan + i,6] = panels[i,3]
                
                # Bottom Right (swapping first and second, and the y-value becomes negative)
                panels[numPanelsSpan + i,7] = panels[i,10]
                panels[numPanelsSpan + i,8] = -panels[i,11]
                panels[numPanelsSpan + i,9] = panels[i,12]
                
                # Bottom Left (swapping first and second, and the y-value becomes negative)
                panels[numPanelsSpan + i,10] = panels[i,7]
                panels[numPanelsSpan + i,11] = -panels[i,8]
                panels[numPanelsSpan + i,12] = panels[i,9]

                # Setting a flag to tell us whether this is a trailing edge panel
                if j == numPanelsChord
                    panels[numPanelsSpan + i,13] = 1
                end
                
            end
            
            # Ordering the panels from most negative Y to most positive Y
            for i = 1:numPanelsSpan

                # Putting the negative Y panels at the front of the array
                orderedPanels[i,:] = panels[end + 1 - i,:]

            end

            for i = 1:numPanelsSpan

                orderedPanels[numPanelsSpan + i,:] = panels[i,:]

            end

            finalPanels = cat(dims=1,finalPanels,orderedPanels)

        end

        return finalPanels
        
    end

    function defineFreestream(panels)
    
        # This function can be changed to determine what goes inside of U
        # Please don't delete any possiblities, just comment out the ones you don't want at the moment
        
        U = 1*ones(length(panels[:,1]));
        
    end

    function createAIC(panels,latticeType)
    
        coordinates = definePoints(panels)
        Xm = coordinates[:,1]
        Ym = coordinates[:,2]
        Zm = coordinates[:,3]
        X1n = coordinates[:,4]
        Y1n = coordinates[:,5]
        Z1n = coordinates[:,6]
        X2n = coordinates[:,7]
        Y2n = coordinates[:,8]
        Z2n = coordinates[:,9]
        unitNormals = zeros(length(panels[:,1]),3)
    
    
        AIC = zeros(length(panels[:,1]),length(panels[:,1])) # Initializing our Aerodynamic Influence Coefficients Matrix
    
        if latticeType == "Horseshoe Lattice"
            
            for i = 1:length(panels[:,1]) # Control Points
            
                # Creating two vectors that are known to lie in the plane of the panel
                panelVector1 = [Xm[i]-X1n[i],Ym[i]-Y1n[i],Zm[i]-Z1n[i]];
                panelVector2 = [Xm[i]-X2n[i],Ym[i]-Y2n[i],Zm[i]-Z2n[i]];
                
                # Crossing those two vectors to find a vector that is known to be normal to the panel
                # The first term is panelVector2 because that is what gives a vector in the positive direction
                normalDirection = cross(panelVector2,panelVector1);
                
                # Normalizing that vector to create the unit normal
                unitNormals[i,:] = normalDirection./norm(normalDirection);
                unitNormals[i,1] = -unitNormals[i,1] # FIXME: I shouldn't have to do this manually to get the x-component in the right direction
                #println(unitNormals[i,:])
    
                for j = 1:length(panels[:,1]) # Influence of each panel on control point "i"
    
                    r1 = [Xm[i],Ym[i],Zm[i]] - [X1n[j],Y1n[j],Z1n[j]];
                    r2 = [Xm[i],Ym[i],Zm[i]] - [X2n[j],Y2n[j],Z2n[j]];
    
                    AIC[i,j] = dot(Velocity(r1,r2,1,"Horseshoe"),unitNormals[i,:]);
            
                end
            
            end
    
        end
    
        if latticeType == "Ring Lattice"
    
            for i = 1:length(panels[:,1]) # Control Points
    
                # We need to find the center of each panel and make that the control point
                X = (panels[i,1] + panels[i,4] + panels[i,7] + panels[i,10]) / 4
                Y = (panels[i,2] + panels[i,5] + panels[i,8] + panels[i,11]) / 4
                Z = (panels[i,3] + panels[i,6] + panels[i,9] + panels[i,12]) / 4
            
                # Creating two vectors that are known to lie in the plane of the panel
                panelVector1 = [Xm[i]-X1n[i],Ym[i]-Y1n[i],Zm[i]-Z1n[i]];
                panelVector2 = [Xm[i]-X2n[i],Ym[i]-Y2n[i],Zm[i]-Z2n[i]];
                
                # Crossing those two vectors to find a vector that is known to be normal to the panel
                # The first term is panelVector2 because that is what gives a vector in the positive direction
                normalDirection = cross(panelVector2,panelVector1);
                
                # Normalizing that vector to create the unit normal
                unitNormals[i,:] = normalDirection./norm(normalDirection);
                unitNormals[i,1] = -unitNormals[i,1] # FIXME: I shouldn't have to do this manually to get the x-component in the right direction
                #println(unitNormals[i,:])
    
                for j = 1:length(panels[:,1]) # Influence of each panel on control point "i"
    
                    # Front vortex
                    r1 = [X,Y,Z] - [panels[j,1],panels[j,2],panels[j,3]];
                    r2 = [X,Y,Z] - [panels[j,4],panels[j,5],panels[j,6]];
                    velocity = Velocity(r1,r2,1,"Bound")
    
                    # Right vortex
                    r1 = [X,Y,Z] - [panels[j,4],panels[j,5],panels[j,6]];
                    r2 = [X,Y,Z] - [panels[j,7],panels[j,8],panels[j,9]];
                    velocity = velocity .+ Velocity(r1,r2,1,"Bound")
    
                    # Rear vortex
                    r1 = [X,Y,Z] - [panels[j,7],panels[j,8],panels[j,9]];
                    r2 = [X,Y,Z] - [panels[j,10],panels[j,11],panels[j,12]];
                    velocity = velocity .+ Velocity(r1,r2,1,"Bound")
    
                    # Left vortex
                    r1 = [X,Y,Z] - [panels[j,10],panels[j,11],panels[j,12]];
                    r2 = [X,Y,Z] - [panels[j,1],panels[j,2],panels[j,3]];
                    velocity = velocity .+ Velocity(r1,r2,1,"Bound")
    
                    # Adding a trailing horseshoe
                    if panels[j,13] == 1 # check the trailing edge flag
    
                        # clear the current velocity, we want a horseshoe vortex, not a ring
                        velocity = [0,0,0]
    
                        # Use the trailing edge of the panel to define a horseshoe vortex that trails off the wing
                        r1 = [Xm[i],Ym[i],Zm[i]] - [X1n[j],Y1n[j],Z1n[j]];
                        r2 = [Xm[i],Ym[i],Zm[i]] - [X2n[j],Y2n[j],Z2n[j]];
                        velocity = velocity .+ Velocity(r1,r2,1,"Horseshoe")
    
                    end
    
                    AIC[i,j] = dot(velocity,unitNormals[i,:]);
    
                end
            
            end
    
        end
        
        # In the AIC, each row represents a control point and each 
        # column is the influence of a different panel on that 
        # control point
    
        return AIC, unitNormals
        
    end

    function calculateLift(density,freestream,panels,GammaValues)
        # Calculates lift, assuming panels are parallelograms
        # This section is already generalized to 3 dimensions because the planform area is defined as the area
        # projected into the x-y plane as viewed from above the wing
        Lift = 0; # Initial lift assumed to be zero, then we'll add the differential lift given by each panel
        CL = 0; # Initial lift coefficient assumed to be zero. We'll fix that later
        cl = zeros(length(panels[:,1])); # Will become our distribution of lift coefficient per panel, ClSpanLocations = zeros(length(panels[:,1])); # Will become our panel center locations
        clSpanLocations = zeros(length(panels[:,1]))
        Area = 0; # Initial planform area assumed to be zero. We'll fix that later

        for i = 1:length(panels[:,1])
        
            deltaY = abs(panels[i,2] - panels[i,5])
            
            Lift = Lift + density*norm(freestream[i,:])*GammaValues[i]*deltaY;
            
        end
        
        dynamicPressure = 0.5*density*norm(freestream[1,:])^2 # FIXME: Make more general. It only uses one panel now
        # We had to take the average freestream velocity so that dynamic pressure would be a scalar
        
        # Finding the total planform area
        for i = 1:length(panels[:,1])
            
            # The area of a parallelogram is equal to its base times its height
            # We will take each deltaY to be our height and use deltaX for our base
            
            deltaY = abs(panels[i,2] - panels[i,5])
            deltaX = abs(panels[i,1] - panels[i,10])
            panelFrontLength = sqrt((panels[i,1] - panels[i,4])^2 + (panels[i,2] - panels[i,5])^2)

            # Find all four sides. We will define them as vectors traveling clockwise around the panel when viewed from above.
            a = [panels[i,4] - panels[i,1], panels[i,5] - panels[i,2], panels[i,6] - panels[i,3]]
            b = [panels[i,7] - panels[i,4], panels[i,8] - panels[i,5], panels[i,9] - panels[i,6]]
            c = [panels[i,10] - panels[i,7], panels[i,11] - panels[i,8], panels[i,12] - panels[i,9]]
            d = [panels[i,1] - panels[i,10], panels[i,2] - panels[i,11], panels[i,3] - panels[i,12]]

            # Find the angles between opposite corners of the quadrilateral
            theta_a_b = acos(dot(a,b) / (norm(a)*norm(b)))
            theta_c_d = acos(dot(c,d) / (norm(c)*norm(d)))

            incrementalArea = 0.5 * norm(a) * norm(b) * sin(theta_a_b) + 0.5 * norm(c) * norm(d) * sin(theta_c_d)
            
            Area = Area + incrementalArea
            
            # See eqn 7.51 in Bertin's book
            cl[i] = (density.*freestream[i]*GammaValues[i]*deltaY)/(dynamicPressure * incrementalArea) # if unit span is only in y
            #Cl[i] = (freestream[i]*GammaValues[i]*deltaY)/(dynamicPressure * deltaX * panelFrontLength) # if unit span is along the leading edge
            clSpanLocations[i] = panels[i,2] + deltaY/2

        end
        
        # Finding the Total Lift Coefficient. We could have added up the Cl values but I'd rather get it in one
        # simple formula with less rounding error
        CL = Lift/(dynamicPressure * Area)
        
        return CL, cl, clSpanLocations;
    end

    function calculateInducedVelocity(panels,GammaValues,location)
    
        boundVortexCenters = zeros(length(panels[:,1]),3); # Will store the center of each bound vortex
        inducedVelocity = zeros(length(panels[:,1])); # Will become our induced velocity at each panel control point
    
        coordinates = definePoints(panels);
        Xm = coordinates[:,1]
        Ym = coordinates[:,2]
        Zm = coordinates[:,3]
        X1n = coordinates[:,4]
        Y1n = coordinates[:,5]
        Z1n = coordinates[:,6]
        X2n = coordinates[:,7]
        Y2n = coordinates[:,8]
        Z2n = coordinates[:,9]
        controlPoints = coordinates[:,1:3];
    
        # Calclulating the center of each bound vortex
        for i = 1:length(panels[:,1])
    
            # deltaX = abs(panels[i,4] - panels[i,1])
            # deltaY = abs(panels[i,5] - panels[i,2])
            # deltaZ = abs(panels[i,6] - panels[i,3])

            deltaX = X2n[i] - X1n[i]
            deltaY = Y2n[i] - Y1n[i]
            deltaZ = Z2n[i] - Z1n[i]
    
            boundVortexCenters[i,1:3] = [X1n[i] + deltaX/2, Y1n[i] + deltaY/2, Z1n[i] + deltaZ/2]
    
        end

        # figure()
        # plotPanels(panels)
        # scatter3D(boundVortexCenters[:,1],boundVortexCenters[:,2],boundVortexCenters[:,3],color = "black")
    
        ###############################
        # At the center of the bound vortex
        if location == "quarter chord"
            # Calculating the induced velocities
            for i = 1:length(panels[:,1])
                currentPoint = boundVortexCenters[i,1:3]
    
                for j = 1:length(panels[:,1]) # Each horseshoe vortex associated with a panel
                
                    r1 = currentPoint - [X1n[j],Y1n[j],Z1n[j]];
                    r2 = currentPoint - [X2n[j],Y2n[j],Z2n[j]];
    
                    if j != i
                        inducedVelocity[i] = inducedVelocity[i] + Velocity(r1,r2,GammaValues[j],"Horseshoe")[3]; # Only uses the z-component of the induced velocity
                    elseif j == i # ignore the bound vortex
                        LeftFilament =  Velocity(r1,r2,GammaValues[j],"Semi-Infinite Left")[3]; # Semi-infinite filament starting at  the lower y-value of the bound vortex
                        RightFilament = Velocity(r1,r2,GammaValues[j],"Semi-Infinite Right")[3]; # Semi-infinite filament starting at  the upper y-value of the bound vortex
                        inducedVelocity[i] = inducedVelocity[i] + LeftFilament + RightFilament
                    end
    
                end
    
            end
        end
    
    
        # ###############################
    
        # ###############################
        # # At the control point at 3/4 chord
    
        # if location == "three quarter chord"
        #     # Calclulating the induced velocities
        #     for i = 1:length(controlPoints[:,1]) # Number of panels
            
        #     currentPoint = [Xm[i],Ym[i],Zm[i]]
        #     # println("Current Point: ",currentPoint)
        #         for j = 1:length(panels[:,1]) # Each horseshoe vortex associated with a panel
                    
        #             r1 = currentPoint - [X1n[j],Y1n[j],Zm[j]];
        #             r2 = currentPoint - [X2n[j],Y2n[j],Zm[j]];
                
        #             inducedVelocity[i] = inducedVelocity[i] + Velocity(r1,r2)[3] * GammaValues[j] / (4*pi); # Only uses the z-component of the induced velocity
                    
        #         end
        #     end 
        # end
    
        # ###############################
    
        # ###############################
        # # Center of the leading edge
    
        # # Calclulating the induced velocities
        # # for i = 1:length(controlPoints[:,1]) # Number of panels
        
        # #     currentPoint = [(panels[i,1] + panels[i,4])/2,(panels[i,2] + panels[i,5])/2,(panels[i,3] + panels[i,6])/2]
        # #     # println("Current Point: ",currentPoint)
        # #     for j = 1:length(panels[:,1]) # Each horseshoe vortex associated with a panel
                
        # #         r1 = currentPoint - [X1n[j],Y1n[j],Zm[j]];
        # #         r2 = currentPoint - [X2n[j],Y2n[j],Zm[j]];
    
        # #         inducedVelocity[i] = inducedVelocity[i] + Velocity(r1,r2)[3] * GammaValues[i] / (4*pi); # Only uses the z-component of the induced velocity
                
        # #     end
        # # end
    
        # ###############################
    
        return inducedVelocity
    
    end

    function calculateInducedDrag(density,freestream,panels,GammaValues,CLn)
        
        # Calculates induced drag, assuming panels are parallelograms
        # This section is already generalized to 3 dimensions because the planform area is defined as the area
        # projected into the x-y plane as viewed from above the wing   

        CDi_near = 0; # Initial induced drag coefficient assumed to be zero. We'll fix that later
        Cd = zeros(length(panels[:,1])); # Will become our distribution of induced drag coefficient per panel
        CdnSpanLocation = zeros(length(panels[:,1])); # Will become our panel center locations
        Area = 0; # Initial planform area assumed to be zero. We'll fix that later
        inducedVelocity = zeros(length(panels[:,1])); # Will become our induced velocity at each panel control point
        inducedAlpha = zeros(length(panels[:,1])); # Will store induced Angle of Attack Information
        boundVortexCenters = zeros(length(panels[:,1]),3); # Will store the center of each bound vortex
        dynamicPressure = 0.5*density*norm(freestream[1,:])^2;
        
        ##############################
        # Center of the bound vortex

        inducedVelocity = calculateInducedVelocity(panels,GammaValues,"quarter chord") 

        #integrate over the lifting line
        for i = 1:length(inducedVelocity)

            deltaY = abs(panels[i,2] - panels[i,5])
            deltaX = abs(panels[i,1] - panels[i,10])

            Cd[i] = -(density * inducedVelocity[i] * GammaValues[i] * deltaY) ./ (dynamicPressure * deltaX * deltaY)

            CdnSpanLocation[i] = panels[i,2] + deltaY/2

            CDi_near = CDi_near + Cd[i]*deltaY # multiply by deltaY because Cd[i] is the induced drag per unit span
            # Does this mean that it should be deltaY/span?
            # See Bertin's Book, Eqn 5.23 on page 218

        end
        
        return CDi_near, Cd, CdnSpanLocation, inducedVelocity

    end

    function calculateDrag(density,freestream,panels,GammaValues)
    
        d = zeros(length(panels[:,1]),1)
        Cd = zeros(length(panels[:,1]),1)
        CdSpanLocations = zeros(length(panels[:,1]),1)
        Drag = 0
        Area = 0
        dynamicPressure = 0.5*density*norm(freestream[1,:])^2
        
        for i = 1:length(panels[:,1]) # Iterating through each panel
            
            yi = panels[i,2] # Front "left" corner from leading edge
            yiPlusOne = panels[i,5] # Front "right" corner from leading edge
                
            zi = panels[i,3] # Front "left" corner from leading edge
            ziPlusOne = panels[i,6] # Front "right" corner from leading edge
            
            yiBar = 0.5 * (yiPlusOne + yi)
            ziBar = 0.5 * (ziPlusOne + zi)
            
            jSummation = 0
            
            for j = 1:length(panels[:,1]) # Iterating through the effects of each other panel on the current panel
                
                yj = panels[j,2] # Front "left" corner from leading edge
                yjPlusOne = panels[j,5] # Front "right" corner from leading edge
                
                zj = panels[j,3] # Front "left" corner from leading edge
                zjPlusOne = panels[j,6] # Front "right" corner from leading edge
                
                kij = (ziBar - zj) * (ziPlusOne - zi) + (yiBar - yj) * (yiPlusOne - yi)
                kij = kij/((yiBar - yj)^2 + (ziBar - zj)^2)
                
                kijPlusOne = (ziBar - zjPlusOne) * (ziPlusOne - zi) + (yiBar - yjPlusOne) * (yiPlusOne - yi)
                kijPlusOne = kijPlusOne/((yiBar - yjPlusOne)^2 + (ziBar - zjPlusOne)^2)
                
                jSummation = jSummation + GammaValues[j] * (kij - kijPlusOne)
    
            end
    
            iSummation = GammaValues[i] * jSummation # Using the second summation (the one over "i")
    
            d[i] = iSummation * density / (4*pi) # distribution of drag
            CdSpanLocations[i] = yiBar # for plotting with Cd
            
            Drag = Drag + d[i] # total drag
    
        end
        
        for i = 1:length(panels[:,1])
            
            # The area of a parallelogram is equal to its base times its height
            # We will take each deltaY to be our height and use deltaX for our base
            
            deltaY = abs(panels[i,2] - panels[i,5])
            deltaX = abs(panels[i,1] - panels[i,10])
            
            incrementalArea = deltaX * deltaY
            
            Area = Area + incrementalArea
            Cd[i] = d[i]/(dynamicPressure * deltaX * deltaY)
            
        end
    
        CDi_far = Drag/(dynamicPressure * Area)
        
        return CDi_far, Cd, CdSpanLocations
        
    end

    function plotPanels(panels)
    
        # Initializing the array
        pointsToPlot = zeros(4*length(panels[:,1]),3)
        
        for i = 1:length(panels[:,1])
            
            pointsToPlot[4*(i-1) + 1,1] = panels[i,1]
            pointsToPlot[4*(i-1) + 1,2] = panels[i,2]
            pointsToPlot[4*(i-1) + 1,3] = panels[i,3]
            pointsToPlot[4*(i-1) + 2,1] = panels[i,4]
            pointsToPlot[4*(i-1) + 2,2] = panels[i,5]
            pointsToPlot[4*(i-1) + 2,3] = panels[i,6]
            pointsToPlot[4*(i-1) + 3,1] = panels[i,7]
            pointsToPlot[4*(i-1) + 3,2] = panels[i,8]
            pointsToPlot[4*(i-1) + 3,3] = panels[i,9]
            pointsToPlot[4*(i-1) + 4,1] = panels[i,10]
            pointsToPlot[4*(i-1) + 4,2] = panels[i,11]
            pointsToPlot[4*(i-1) + 4,3] = panels[i,12]
            
        end
    
        for i = 1:Int(length(pointsToPlot[:,1])/4)
            
            firstIndex = Int(4*(i-1)+1)
            secondIndex = Int(firstIndex + 3)
            
            PyPlot.surf(pointsToPlot[firstIndex:secondIndex,1],
                 pointsToPlot[firstIndex:secondIndex,2],
                 pointsToPlot[firstIndex:secondIndex,3]) 
            PyPlot.zlim(-.5,.5)
            
        end
        
    end

    function plotLiftDistribution(panels,cl)
    
        # Initializing the array
        pointsToPlot = zeros(4*length(panels[:,1]),3)
        
        for i = 1:length(panels[:,1])
            
            pointsToPlot[4*(i-1) + 1,1] = panels[i,1]
            pointsToPlot[4*(i-1) + 1,2] = panels[i,2]
            pointsToPlot[4*(i-1) + 1,3] = panels[i,3]
            pointsToPlot[4*(i-1) + 2,1] = panels[i,4]
            pointsToPlot[4*(i-1) + 2,2] = panels[i,5]
            pointsToPlot[4*(i-1) + 2,3] = panels[i,6]
            pointsToPlot[4*(i-1) + 3,1] = panels[i,7]
            pointsToPlot[4*(i-1) + 3,2] = panels[i,8]
            pointsToPlot[4*(i-1) + 3,3] = panels[i,9]
            pointsToPlot[4*(i-1) + 4,1] = panels[i,10]
            pointsToPlot[4*(i-1) + 4,2] = panels[i,11]
            pointsToPlot[4*(i-1) + 4,3] = panels[i,12]
            
        end
    
        for i = 1:Int(length(pointsToPlot[:,1])/4)
            
            firstIndex = Int(4*(i-1)+1)
            secondIndex = Int(firstIndex + 3)
            
            # Grayscale colors: color = string(1 - cl[i]/maximum(cl))
    
            PyPlot.surf(pointsToPlot[firstIndex:secondIndex,1],
                 pointsToPlot[firstIndex:secondIndex,2],
                 pointsToPlot[firstIndex:secondIndex,3],color = [cl[i]/maximum(cl),1-cl[i]/maximum(cl),1-cl[i]/maximum(cl)]) 
            PyPlot.zlim(-.5,.5)
            
        end
        
    end

    function plotInducedDragDistribution(panels,cd)
    
        # Initializing the array
        pointsToPlot = zeros(4*length(panels[:,1]),3)
        
        for i = 1:length(panels[:,1])
            
            pointsToPlot[4*(i-1) + 1,1] = panels[i,1]
            pointsToPlot[4*(i-1) + 1,2] = panels[i,2]
            pointsToPlot[4*(i-1) + 1,3] = panels[i,3]
            pointsToPlot[4*(i-1) + 2,1] = panels[i,4]
            pointsToPlot[4*(i-1) + 2,2] = panels[i,5]
            pointsToPlot[4*(i-1) + 2,3] = panels[i,6]
            pointsToPlot[4*(i-1) + 3,1] = panels[i,7]
            pointsToPlot[4*(i-1) + 3,2] = panels[i,8]
            pointsToPlot[4*(i-1) + 3,3] = panels[i,9]
            pointsToPlot[4*(i-1) + 4,1] = panels[i,10]
            pointsToPlot[4*(i-1) + 4,2] = panels[i,11]
            pointsToPlot[4*(i-1) + 4,3] = panels[i,12]
            
        end
    
        for i = 1:Int(length(pointsToPlot[:,1])/4)
            
            firstIndex = Int(4*(i-1)+1)
            secondIndex = Int(firstIndex + 3)
            
            # Grayscale colors: color = string(1 - cl[i]/maximum(cl))
            println(abs(cd[i]/maximum(cd)))
    
            PyPlot.surf(pointsToPlot[firstIndex:secondIndex,1],
                 pointsToPlot[firstIndex:secondIndex,2],
                 pointsToPlot[firstIndex:secondIndex,3],color = abs(cd[i]/maximum(cd))) 
            PyPlot.zlim(-.5,.5)
            
        end
        
    end

end