% Purpose: This function takes a static fire test CSV file and completes a
% full basic analysis of the data from the test.
%
% Input:
%   csvFile = a CSV file that contains the test information and is
%             formatted the same way as "Static Fire Test.csv" in the 
%             project on GitLab.
%
% Output:
%   A PDF file that contains the data segarated as:
%       - Array-wise: OASPL, Derivative Skewness
%       - Each Station: Waveform, OASPL, Derivative Skewness, Spectrum
%       - Weather Report: Wind Magnitude, Wind Magnitude + Direction
%


function StaticFire(csvFile)

    clearvars -except csvFile

    % You guys will need to change the following two paths to match what's
    % on your computers. It's a temporary fix. We'll find a way to automate
    % the path.
    
    % This is the path to the "General Signal Processing Library" that you
    % can clone from GitLab under BYU Acoustics > General Acoustics
    % Research > General Signal Processing Library
    % addpath("/Users/markanderson/Box/Acoustics-MCA/Code/General Acoustics Research/General Signal Processing Library")
    
    % This is the path to the data on your computer. I've shared the sample
    % data with you on Box
    pathToData = '/Users/markanderson/Box/Acoustics-MCA/Other/Projects/Sample Static Fire Data/GEM-63 October 2019';

    % From here on down, it should run without you needing to change any
    % paths.
    disp('Extracting CSV File Contents')
    [microphoneInfo, sourceInfo, analysisOptions,arrayResults, stationResults] = ...
        parseCSV(csvFile);
    
    % Extract all of the waveforms
    disp('Extracting the Data')
    [longWaveforms,shortWaveforms] = ... 
        getWaveform(microphoneInfo,analysisOptions,pathToData);
    
    % Calculate the OASPL values
    disp('Calculating OASPL Values')
    [everyOASPL,rollingOASPL,stationOASPL] = ...
        getOASPL_StaticFire(longWaveforms,shortWaveforms,microphoneInfo);
    
    % Calculate the Derivative Skewness values
    disp('Calculating Derivative Skewness Values')
    [everyDerivativeSkewness,rollingDerivativeSkewness,...
          stationDerivativeSkewness] = ...
          getDerivativeSkewness_StaticFire(longWaveforms,shortWaveforms,microphoneInfo);
    
    %-- Creating the Report --%
    % Making the Title Page
    titleText = ["Analysis of the GEM-63 Static Fire Test",...
                 "",...
                 "Location: Promontory, Utah",...
                 "",...
                 "Date: Septemeber, 2019"];
    titlePage(titleText)
    
    titleText = ["Array-wise Data"];
    titlePage(titleText)
    
    % Plotting along the array and saving it to a PS file
    plotArray(microphoneInfo,everyOASPL,stationOASPL,...
         everyDerivativeSkewness,stationDerivativeSkewness,shortWaveforms)  
      
    titleText = ["Data At Each Station"];
    titlePage(titleText) 
     
    % Plotting each station and saving it to a PS file
    plotStations(longWaveforms,shortWaveforms,rollingOASPL,...
         rollingDerivativeSkewness,microphoneInfo)
     
    titleText = ["Weather Report"];
    titlePage(titleText) 
      
    close all
    clear titleText
    
    %disp('Saving the workspace')
    %save('Output/Static Fire Data.mat','-v7.3')

end % StaticFire