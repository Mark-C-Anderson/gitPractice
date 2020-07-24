module makeAirfoil

    using Pkg
    Pkg.add(PackageSpec(url="https://github.com/byuflowlab/Xfoil.jl"))
    Pkg.add("CSV")
    Pkg.add("Tables")
    import Xfoil
    using CSV
    using Tables
    using AbstractAlgebra

    function makeXarray(chord,dphi)

        phi = 0:dphi:pi
        x = zeros(length(phi),1)
        for i = 1:length(phi)
            x[i] = chord/2 * (1-cos(phi[i]))
        end

        return x

    end # makeXarray

    function naca(c,p,t,dxi)
        
        c = c/100;
        p = p/10;
        t = t/100;
 
        x = makeXarray(1,dphi)
        array = zeros(length(x),3);
        coordinates = zeros(2*length(x),2);

        for i = 1:1:length(x);
            array[i,1] = x[i];
            array[i,2] = camber(x[i],p,c) + thickness(t,x[i])/2;
            array[i,3] = camber(x[i],p,c) - thickness(t,x[i])/2;
        end

        xpositions = Vector(array[:,1]);
        zu = Vector(array[:,2]);
        zl = Vector(array[:,3]);

        coordinates[1:length(x),1] = xpositions;
        coordinates[length(x) + 1:2*length(x),1] = reverse(xpositions);
        coordinates[1:length(x),2] = zu;
        coordinates[length(x) + 1:2*length(x),2] = reverse(zl);

        # XFOIL wants the trailing edge first, then wrap around to the leading edge
        # and then back to the trailing edge
        # Currently, we have it:
        # leading edge -> trailing edge over the top
        # leading edge -> trailing edge over the bottom
        # So all we need to do is reverse the top, or the first half of the coordinates

        coordinates[1:length(x),1] = reverse(coordinates[1:length(x),1])
        coordinates[1:length(x),2] = reverse(coordinates[1:length(x),2])

        coordinates[length(x) + 1:2*length(x),1] = reverse(coordinates[length(x) + 1:2*length(x),1])
        coordinates[length(x) + 1:2*length(x),2] = reverse(coordinates[length(x) + 1:2*length(x),2])

        runXFOIL(a,b,c)

        return coordinates

    end # naca

    function camber(x,p,c)

        # camber(x,p,c)

        # Determines the proper amount of camber.
        
        # x = position along the chord
        # p = position of maximum camber
        # c = value of maximum camber

        if x <= p
            z = c * (2*p*x - x^2)/p^2;
        else
            z = c * (1 - 2*p + 2*p*x -x^2)/(1-p)^2;
        end
    
        # BUT... if p = 0
        # if p == 0 && x == 0
        #     z = 0
        # elseif p == 0 && x == 1
        #     z = 0
        # end
    
        if p == 0
            z = 0
        end
    
        return z
    
    end

    function thickness(m,x)

        # thickness(maximum thickness,position)

        # Outputs the thickness of a NACA 4-series airfoil.
        
        # This function the maximum thickness (as a percentage of the chord) 
        # and the position along the chord to generate the thickness along
        # the chord.

        t = 10 * m * (0.2969*sqrt(x) - 0.1260*x - 0.3537*x^2 + 0.2843*x^3 - 0.1015x^4);
        
    end

    function tabulateData(airfoil,angleRange,reynoldsNumberRange,path,airfoilName)

        totalPath = string(path,"/",airfoilName,".csv")

        println("Calculating data...")

        # Printing out the directory where the data are saved
        println("Data will be stored to:")
        println(totalPath)

        # Creating the .csv file
        touch(totalPath)

        # Performing the analysis
        data = zeros(length(angleRange) + 1, length(reynoldsNumberRange) + 1)
        data[2:end,1] = angleRange
        data[1,2:end] = reynoldsNumberRange
        reynoldsNumber = reynoldsNumberRange[1]

        for i = 1:length(reynoldsNumberRange)
            reynoldsNumber = reynoldsNumberRange[i]
            println("Calculating data for Re: ", reynoldsNumber)
            cl, cdd, cdp, cm, converged = Xfoil.xfoilsweep(airfoil[:,1],airfoil[:,2],angleRange,reynoldsNumber)
            data[2:end,i+1] = cl
        end

        # Writing to the .csv file
        data = Tables.table(data)
        CSV.write(totalPath,data)

    end # tabulateData

end # module
