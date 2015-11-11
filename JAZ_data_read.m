% JAZ_data_read.m
%
% Created 5 Nov 2013 by JRC under MATLAB R2013a; version history subsequent
% to 5 Nov 2015 on GitHub
% Owner/author: James Collins, MIT-WHOI Joint Program, Woods Hole
% Oceanographic Institution, james.r.collins@aya.yale.edu
%
% Purpose: Read in and store to file UV-VIS spectra from the Ocean Optics
% JAZ device
%
% Dependencies and required files:
%
% 1. A deployment log file (example maintained in same GitHub repo where you got
% this from; Jaz_deployment_log_example.xlsx)
% 
% 2. For subsequent analysis, the .csv file JAZ_wavelengths.csv (also
% available in the repo), which contains a list of the wavelengths measured
% by the instrument (list provided specific to the instrument J.R.C. used
% to collect data; verify that the range of your instrument is same)

%% Read in data

% Specify the directory into which you've downloaded the data from the JAZ
% SD card; the folders in this directory representing each data download
% should have the format "Download_YYYYMMDD"; these folders should contain
% the "irrad_" subfolders that you've copied directly off the device SD
% card (without any further modification)

JAZfiles_directory='/Volumes/dept/C-045/PAL 15-16/PAR and UV-VIS Spec Data/UV-VIS data from Jaz spectrometer (Ocean Optics, Inc.)/';

% Query to find only folders that ultimately contain JAZ data
JAZDownloadfolders_only=strcat(JAZfiles_directory,'Download*');

% Execute query to find the folders that have spectral data ("irrad" directories) in them 
JAZDownloadfolders = dir(JAZDownloadfolders_only);

insertrow=1; % Set to 1 to build a new file starting at row 1, or can manually
% change to any row in your JAZdata matrix where you'd like new data
% appended

suspectcount=1; % To keep track of how many suspect records we have

for i=1:length(JAZDownloadfolders)
       thisdir = JAZDownloadfolders(i,1).name; % Get the name of the directory
       thisdirpath = strcat(JAZfiles_directory,thisdir,'/');
       irraddirs_only = strcat(thisdirpath,'irrad*'); % Query to find only the subdirectories
       % that have files from the JAZ in them
       irraddirs = dir(irraddirs_only); % Execute the query using dir()
       for j=1:length(irraddirs)
                  thissubdir = irraddirs(j,1).name; % Get the name of the directory
                  thissubdirpath = strcat(JAZfiles_directory,thisdir,'/',thissubdir,'/');
                  outputfiles_only = strcat(thissubdirpath,'OUTPUTFILE*');
                  % Query to find only JAZ "OUTPUTXXXX" files
                  outputfiles = dir(outputfiles_only); % Execute the query using dir()
                  for k=1:length(outputfiles)
                      thisoutputfile = outputfiles(k,1).name; % Get the name of the file
                      thisoutputfilepath = strcat(JAZfiles_directory,thisdir,...
                          '/',thissubdir,'/',thisoutputfile);
                      % Extract the spectrum timestamp in GMT
                      fileID = fopen(thisoutputfilepath);
                      formatspec_timestamp = '%*s %24c';
                      spectrum_time_raw = textscan(fileID,formatspec_timestamp,1,...
                          'HeaderLines',2);
                      fclose(fileID);
                      % Convert date & time in GMT to MATLAB-speak
                      formatIn = 'ddd mmm dd HH:MM:SS yyyy';
                      spectrum_time = datenum(spectrum_time_raw,formatIn);
                      % Extract the integration time
                      fileID = fopen(thisoutputfilepath);
                      formatspec_inttime = '%*s %*s %*s %f';
                      inttime_raw = textscan(fileID,formatspec_inttime,1,...
                          'HeaderLines',7);
                      fclose(fileID);
                      inttime = cell2mat(inttime_raw);
                      % Collect data
                      fileID = fopen(thisoutputfilepath)
                      formatspec_specdata = '%f %f %f %f';
                      specdata_raw = textscan(fileID,formatspec_specdata,...
                          'HeaderLines',20);
                      fclose(fileID);
                      specdata = cell2mat(specdata_raw);
                      % Do a little test to see if the instrument was
                      % saturated at any wavelength during this
                      % measurement; we can infer this if there are many
                      % repeated high values in col 3 of the JAZ data output
                      % we'll use 5 repetitions as the threshold to flag
                      % the spectrum for further investigation, and write
                      % some information about those suspect spectra to a
                      % matrix "invest_for_possible_saturation"
                      [a,b]=hist(specdata(:,3),unique(specdata(:,3)));
                      if max(a) >= 5
                          invest_for_possible_saturation(suspectcount,1:3) = ...
                              [spectrum_time inttime max(a)];
                          suspectcount=suspectcount+1;
                      end
                      % Now, write data to a matrix; only going to write
                      % the final, processed irradiances from col 4 of the
                      % JAZ output file
                      JAZdata(insertrow,1:2051) = [spectrum_time inttime max(a) specdata(:,4)'];
                      % Update plot as we read in data
                      figure(1);
                      plot(specdata(212:1756,1),specdata(212:1756,4));
                      title(datestr(spectrum_time));
                      xlabel('Wavelength (nm)');
                      ylabel('Irradiance (uW/cm2/nm)');
                      insertrow=insertrow+1;
                  end
       end
end

%% Make some sense of the suspect, i.e., potentially oversaturated readings

figure;

fig=bar(invest_for_possible_saturation(:,1),invest_for_possible_saturation(:,3));
xlabel('Date');
ylabel('Number of wavelengths at which CCD was saturated');
datetick('x');

%% General analysis

%% For time-series data on a given day, e.g., 20 Nov 13

subset_ind=find(JAZdata(:,1)>=datenum(2013,11,20) & JAZdata(:,1)<=datenum(2013,11,21));
specdata_20Nov=JAZdata(subset_ind,:);
times=specdata_20Nov(:,1);

% Read in wavelengths (in nm)

JAZ_wavelengths = csvread('/Users/jrcollins/Dropbox/High-Lat Lipid Peroxidation/Data/JAZ UV-VIS/JAZ_wavelengths.csv');
JAZ_wavelengths = JAZ_wavelengths';

% Define some spectral ranges (in nm)

UVB=[280 315];
UVA=[315 400];

ind_UVB=find(JAZ_wavelengths>=UVB(1) & JAZ_wavelengths<UVB(2));
ind_UVA=find(JAZ_wavelengths>=UVA(1) & JAZ_wavelengths<UVA(2));

UVB_wavelengths=JAZ_wavelengths(ind_UVB);
UVA_wavelengths=JAZ_wavelengths(ind_UVA);

specdata_20Nov_UVB=specdata_20Nov(:,ind_UVB+3);
specdata_20Nov_UVA=specdata_20Nov(:,ind_UVA+3);
% Offset because there are 3 columns of nonspectral data in the
% specdata_20Nov matrix

% Generate integrated spectral dosages (in uW/cm2) for UVA, UVB at each JAZ timepoint

UVB_dose=nan(length(times),1);
UVA_dose=nan(length(times),1);

for i=1:length(times)
    UVB_dose(i,1) = trapz(UVB_wavelengths,specdata_20Nov_UVB(i,:));
    UVA_dose(i,1) = trapz(UVA_wavelengths,specdata_20Nov_UVA(i,:));
end

% Generate cumulative (time-integrated) dosages (in kJ/m^2) at each timepoint

dt=300; % Time interval = 300 s (this should be whatever sampling interval
        % you used when you collected the data (a JAZ setting) 

sampinv=[1:dt:dt*length(times)];

UVB_dose_cum=nan(length(times),1);
UVA_dose_cum=nan(length(times),1);

for i=2:length(times)
    UVB_dose_cum(i,1) = trapz(sampinv(1:i),UVB_dose(1:i))/10^5;
    UVA_dose_cum(i,1) = trapz(sampinv(1:i),UVA_dose(1:i))/10^5;
end

subplot(2,1,1)
plot(times,UVA_dose,'r-',times,UVB_dose,'b-')
subplot(2,1,2)
plot(times,UVA_dose_cum,'r-',times,UVB_dose_cum,'b-')

dosages_20Nov = [m2xdate(times) UVA_dose UVA_dose_cum UVB_dose UVB_dose_cum];
