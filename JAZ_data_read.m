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

%% Some general treatment of wavelengths

% Read in wavelengths (in nm)

JAZ_wavelengths = csvread('/Users/jrcollins/Code/Optics_Photochem/JAZ_wavelengths.csv');
JAZ_wavelengths = JAZ_wavelengths';

% Define some spectral ranges (in nm)

UVB=[290 315];
UVA=[315 400];

ind_UVB=find(JAZ_wavelengths>=UVB(1) & JAZ_wavelengths<UVB(2));
ind_UVA=find(JAZ_wavelengths>=UVA(1) & JAZ_wavelengths<UVA(2));

UVB_wavelengths=JAZ_wavelengths(ind_UVB);
UVA_wavelengths=JAZ_wavelengths(ind_UVA);

%% Read in data

% Specify the directory into which you've downloaded the data from the JAZ
% SD card; the folders in this directory representing each data download
% should have the format "Download_YYYYMMDD"; these folders should contain
% the "irrad_" subfolders that you've copied directly off the device SD
% card (without any further modification)

% -------- For entire data sets

JAZfiles_directory='/Users/jrcollins/Code/LipidPhotoOxBox/data/raw/JAZ_UV_VIS/raw_spectra/PAL1314/';

% Query to find only folders that ultimately contain JAZ data
JAZDownloadfolders_only=strcat(JAZfiles_directory,'Download*');

% -------- end section for entire data sets

% % -------- For individual wavelength profiles
% 
% JAZfiles_directory='/Users/jrcollins/Code/LipidPhotoOxBox/data/raw/JAZ_UV-VIS/raw_spectra/PAL1516/';
% 
% % Query to find only folders that ultimately contain JAZ data
% JAZDownloadfolders_only=strcat(JAZfiles_directory,'Station B - 20151215*');
% 
% % -------- end section for individual wavelength profiles

% Execute query to find the folders that have spectral data ("irrad" directories) in them 
JAZDownloadfolders = dir(JAZDownloadfolders_only);

insertrow=1; % Set to 1 to build a new file starting at row 1, or can manually
% change to any row in your JAZdata matrix where you'd like new data
% appended

% Counters to keep track of how many suspect records we have
suspectcount=1; 
suspectcount_UVB=1;

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
                      
                      % Do some tests to see if the instrument was
                      % saturated at any wavelength during this
                      % measurement, or saturated in the UVB range;
                      % we can infer this if there are many
                      % repeated values of large numbers
                      % in col 3 of the JAZ data output
                      
                      % We'll use 3 repetitions as the threshold to flag
                      % the spectrum for further investigation, and write
                      % some information about those suspect spectra to
                      % matrices "invest_for_possible_saturation" and 
                      % "invest_for_possible_saturation_UVB"
                      
                      % full spectrum
                      
                      [a_full,b_full]=hist(specdata(:,3),unique(specdata(:,3)));
                      if any(b_full(a_full>=3)==max(specdata(:,3)))
                          invest_for_possible_saturation(suspectcount,1:3) = ...
                              [spectrum_time inttime max(a_full)];
                          suspectcount=suspectcount+1;
                      end
                      
                      % UVB only
                      
                      [a_UVB,b_UVB]=hist(specdata(ind_UVB,3),unique(specdata(ind_UVB,3)));
                      if any(b_UVB(a_UVB>=3)==max(specdata(:,3)))
                          invest_for_possible_saturation_UVB(suspectcount_UVB,1:3) = ...
                              [spectrum_time inttime max(a_UVB)];
                          suspectcount_UVB=suspectcount_UVB+1;
                      end
                      
                      % create an alternate data string 'specdata_QA' in
                      % which the bad data points that reflect CCD
                      % saturation are replaced with 'NaN' 
                      
                      specdata_QA = specdata;
                      
                      for z=1:length(a_full)
                          if (a_full(z)>=3 && b_full(z)==max(specdata(:,3)))
                              specdata_QA(specdata_QA(:,3)==b_full(z),4)=NaN;
                          end
                      end
                      
                      % Now, write data to 2 matrices; only going to write
                      % the final, processed irradiances from col 4 of the
                      % JAZ output file
                      
                      % Without file/scan ID name
                      
                      JAZdata(insertrow,1:2052) = [spectrum_time inttime max(a_full) max(a_UVB) specdata(:,4)'];
                      JAZdata_QA(insertrow,1:2052) = [spectrum_time inttime max(a_full) max(a_UVB) specdata_QA(:,4)'];

%                       % With file/scan ID name
%                       
%                       % Strip out just the scan file ID number out of the
%                       % file name
%                       thisoutputfile_ID = textscan(thisoutputfile,'%10s %4c');
%                       thisoutputfile_ID = str2double(thisoutputfile_ID(2));
% 
%                       JAZdata(insertrow,1:2053) = [thisoutputfile_ID spectrum_time inttime max(a_full) max(a_UVB) specdata(:,4)'];
%                       JAZdata_QA(insertrow,1:2053) = [thisoutputfile_ID spectrum_time inttime max(a_full) max(a_UVB) specdata_QA(:,4)'];
                      
%                       % Update plot as we read in data
%                       figure(1);
%                       plot(specdata(212:1756,1),specdata(212:1756,4));
%                       title(datestr(spectrum_time));
%                       xlabel('Wavelength (nm)');
%                       ylabel('Irradiance (uW/cm2/nm)');

                      insertrow=insertrow+1;
                  end
       end
end

%% save raw data
% full spectrum, not QA'd, all dates, spectra captured both in situ at 0.6
% m water depth and in open air

save('JAZ_UV-VIS_full_spectra_all_PAL1314_uW_cm2.mat','JAZdata')

% %% save data from a single profile
% % full spectra, QA'd
% 
% save('JAZ_UV-VIS_full_spectra_AH_profile_20151215_Stn_B_PAL1516_uW_cm2.mat','JAZdata_QA')
% 
% % export data to .csv; have to use dlmwrite instead of csvwrite to get
% % proper precision
% 
% dlmwrite('JAZ_UV-VIS_full_spectra_AH_profile_20151215_Stn_B_PAL1516_uW_cm2.csv', JAZdata_QA, 'delimiter', ',', 'precision', 15); 

%% Make some sense of the suspect, i.e., potentially oversaturated readings

figure;

fig=bar(invest_for_possible_saturation(:,1),invest_for_possible_saturation(:,3));
xlabel('Date');
ylabel('Number of wavelengths at which CCD was saturated, full spectrum');
datetick('x');

figure;

fig=bar(invest_for_possible_saturation_UVB(:,1),invest_for_possible_saturation_UVB(:,3));
xlabel('Date');
ylabel('Number of wavelengths at which CCD was saturated, UVB wavelengths only');
datetick('x');

% for the PAL1314 data, suspectcount_UVB should be = 1 (and 
% 'invest_for_possible_saturation_UVB' nonexistent at this point in the 
% script), indicating that there were no spectra for which the CCD was 
% saturated in the UVB range

%% General analysis

%% Extract and save "good" in situ spectra

% eliminates those taken when device was deployed out of water, or with
% other issues

Insitu_spectra_PAL1314_uW_cm2_all = JAZdata_QA;

% eliminate readings with missing timestamp

Insitu_spectra_PAL1314_uW_cm2_all = Insitu_spectra_PAL1314_uW_cm2_all(Insitu_spectra_PAL1314_uW_cm2_all(:,1)>0,:);

% extract only "good" data based on JAZ deployment log

Insitu_spectra_PAL1314_uW_cm2_good = ...
    Insitu_spectra_PAL1314_uW_cm2_all(...
    (Insitu_spectra_PAL1314_uW_cm2_all(:,1)<datenum(2013,11,13,9,0,0)) | ...
    (Insitu_spectra_PAL1314_uW_cm2_all(:,1)>datenum(2013,11,14,9,0,0) & ...
    Insitu_spectra_PAL1314_uW_cm2_all(:,1)<datenum(2013,12,3,13,0,0)) | ...
    (Insitu_spectra_PAL1314_uW_cm2_all(:,1)>datenum(2013,12,14,9,30,0) & ...
    Insitu_spectra_PAL1314_uW_cm2_all(:,1)<datenum(2013,12,21,11,0,0)),:);

% % extract data from open-air deployment from 11-13 Nov 13
% 
% Open_air_spectra_PAL1314_11_13Nov13_uW_cm2 = Insitu_spectra_PAL1314_uW_cm2_all(...
%     (Insitu_spectra_PAL1314_uW_cm2_all(:,1)>datenum(2013,11,11,12,0,0) &...
%     Insitu_spectra_PAL1314_uW_cm2_all(:,1)<datenum(2013,11,13,12,0,0)),:);

% quick plot to verify

figure;

plot(Insitu_spectra_PAL1314_uW_cm2_good(:,1),Insitu_spectra_PAL1314_uW_cm2_good(:,30))
datetick('x')     

% save .mat file

save('JAZ_UV-VIS_full_spectra_0.6m_subsurface_PAL1314_uW_cm2_QA.mat','Insitu_spectra_PAL1314_uW_cm2_good')

% export full spectrum data to .csv; have to use dlmwrite instead of csvwrite to get
% proper precision

dlmwrite('JAZ_UV-VIS_full_spectra_0.6m_subsurface_PAL1314_uW_cm2_QA.csv', Insitu_spectra_PAL1314_uW_cm2_good, 'delimiter', ',', 'precision', 15); 

% % export open-air 11-13 Nov 13 data to .csv; have to use dlmwrite instead of csvwrite to get
% % proper precision
% 
% dlmwrite('JAZ_UV-VIS_open_air_full_spectra_PAL1314_11_13Nov13_uW_cm2.csv', Open_air_spectra_PAL1314_11_13Nov13_uW_cm2, 'delimiter', ',', 'precision', 15); 

%% Create subset of just UVB band

Insitu_UVB_spectra_PAL1314_uW_cm2_good=Insitu_spectra_PAL1314_uW_cm2_good(:,[1:4,ind_UVB+4]);

% save .mat file

save('JAZ_UV-VIS_UVB_wavelengths_0.6m_subsurface_PAL1314_uW_cm2_QA.mat','Insitu_UVB_spectra_PAL1314_uW_cm2_good')

% export UVB data to .csv; have to use dlmwrite instead of csvwrite to get
% proper precision

dlmwrite('JAZ_UV-VIS_UVB_wavelengths_0.6m_subsurface_PAL1314_uW_cm2_QA.csv', Insitu_UVB_spectra_PAL1314_uW_cm2_good, 'delimiter', ',', 'precision', 15); 

%% Create subset of just UVA band

Insitu_UVA_spectra_PAL1314_uW_cm2_good=Insitu_spectra_PAL1314_uW_cm2_good(:,[1:4,ind_UVA+4]);

% save .mat file

save('JAZ_UV-VIS_UVA_wavelengths_0.6m_subsurface_PAL1314_uW_cm2_QA.mat','Insitu_UVA_spectra_PAL1314_uW_cm2_good')

% export UVB data to .csv; have to use dlmwrite instead of csvwrite to get
% proper precision

dlmwrite('JAZ_UV-VIS_UVA_wavelengths_0.6m_subsurface_PAL1314_uW_cm2_QA.csv', Insitu_UVA_spectra_PAL1314_uW_cm2_good, 'delimiter', ',', 'precision', 15); 

%% Integrated UVB band calculations

% wavelength-integrated fluxes

UVB_flux_PAL1314_uW_cm2 = nan(size(Insitu_UVB_spectra_PAL1314_uW_cm2_good,1),2);
UVB_flux_PAL1314_uW_cm2(:,1) = Insitu_UVB_spectra_PAL1314_uW_cm2_good(:,1);

for i=1:size(UVB_flux_PAL1314_uW_cm2,1)
    UVB_flux_PAL1314_uW_cm2(i,2) = trapz(UVB_wavelengths,Insitu_UVB_spectra_PAL1314_uW_cm2_good(i,5:end));
end

% eliminate bad data points (those with missing timestamp)

UVB_flux_PAL1314_uW_cm2 = UVB_flux_PAL1314_uW_cm2(UVB_flux_PAL1314_uW_cm2(:,1)>0,:);

% a quick plot

plot(UVB_flux_PAL1314_uW_cm2(:,1),UVB_flux_PAL1314_uW_cm2(:,2))
datetick('x')

% daily doses

% first, for a given day, need to determine whether we have complete 24-hr
% data

% list of all dates for which we have data

UVBdates_all = unique(datetime(year(UVB_flux_PAL1314_uW_cm2(:,1)),...
    month(UVB_flux_PAL1314_uW_cm2(:,1)),...
    day(UVB_flux_PAL1314_uW_cm2(:,1))));

% create date subset based on inspection of plot and info recorded in
% deployment log; eliminates dates for which we have incomplete in-water
% data (either instrument was not deployed at all, or was not deployed at
% 0.6 m water depth)

UVBdates_good = UVBdates_all([2:10,13:15,20:37,41:45]);

% calculate daily integrated fluxes using same convention as NOAA ESRL, see
% http://esrl.noaa.gov/gmd/grad/antuv/docs/netOps/CHAPTER4.PDF, p. 4-25

% preallocate destination matrix

UVB_daily_dose_0_6m_subsurface_PAL1314_kJ_m2 = nan(length(min(UVBdates_all):max(UVBdates_all)),2);
UVB_daily_dose_0_6m_subsurface_PAL1314_kJ_m2(:,1) = datenum(min(UVBdates_all):max(UVBdates_all));

% make calculations with center of each integration period at local noon
% for Palmer, NOAA protocol (link above) defines this as approx. 1600 UTC
%
% note that the JAZ data J.R.C. recorded during the 2013-2014 field season
% at Palmer have timestamps in UTC, *not* local time

% identities

W_per_uW = 1/1000000;
J_per_kJ = 1/1000;
cm2_per_m2 = 10000;

for i=1:size(UVB_daily_dose_0_6m_subsurface_PAL1314_kJ_m2,1)
    if ismember(UVB_daily_dose_0_6m_subsurface_PAL1314_kJ_m2(i,1),datenum(UVBdates_good))
        % define the local noon for this date in MATLAB Julian format
        thisLocalNoon = datenum(UVB_daily_dose_0_6m_subsurface_PAL1314_kJ_m2(i,1)+datenum(0,0,0,16,0,0));
        % select the data within the local noon +/- 12 hr window
        theseUVBData = ...
            UVB_flux_PAL1314_uW_cm2(...
            UVB_flux_PAL1314_uW_cm2(:,1)>(thisLocalNoon-datenum(0,0,0,12,0,0)) & ...
            UVB_flux_PAL1314_uW_cm2(:,1)<(thisLocalNoon+datenum(0,0,0,12,0,0)),:);
        % calculate the time integrated-dose, in kJ/m2
        UVB_daily_dose_0_6m_subsurface_PAL1314_kJ_m2(i,2) = ...
            trapz(theseUVBData(:,1)*24*60*60,theseUVBData(:,2))*...
            W_per_uW*J_per_kJ*cm2_per_m2;
    end
end

% save .mat file

save('Daily_int_UVB_dose_0.6m_subsurface_PAL1314_kJ_m2.mat','UVB_daily_dose_0_6m_subsurface_PAL1314_kJ_m2')

% export UVB doses to .csv; have to use dlmwrite instead of csvwrite to get
% proper precision

dlmwrite('Daily_int_UVB_dose_0.6m_subsurface_PAL1314_kJ_m2.csv', UVB_daily_dose_0_6m_subsurface_PAL1314_kJ_m2, 'delimiter', ',', 'precision', 8); 

%% Integrated UVA band calculations

% wavelength-integrated fluxes

UVA_flux_PAL1314_uW_cm2 = nan(size(Insitu_UVA_spectra_PAL1314_uW_cm2_good,1),2);
UVA_flux_PAL1314_uW_cm2(:,1) = Insitu_UVA_spectra_PAL1314_uW_cm2_good(:,1);

for i=1:size(UVA_flux_PAL1314_uW_cm2,1)
    UVA_flux_PAL1314_uW_cm2(i,2) = trapz(UVA_wavelengths,Insitu_UVA_spectra_PAL1314_uW_cm2_good(i,5:end));
end

% eliminate bad data points (those with missing timestamp)

UVA_flux_PAL1314_uW_cm2 = UVA_flux_PAL1314_uW_cm2(UVA_flux_PAL1314_uW_cm2(:,1)>0,:);

% a quick plot

plot(UVA_flux_PAL1314_uW_cm2(:,1),UVA_flux_PAL1314_uW_cm2(:,2))
datetick('x')

% daily doses

% first, for a given day, need to determine whether we have complete 24-hr
% data

% list of all dates for which we have data

UVAdates_all = unique(datetime(year(UVA_flux_PAL1314_uW_cm2(:,1)),...
    month(UVA_flux_PAL1314_uW_cm2(:,1)),...
    day(UVA_flux_PAL1314_uW_cm2(:,1))));

% create date subset based on inspection of plot and info recorded in
% deployment log; eliminates dates for which we have incomplete in-water
% data (either instrument was not deployed at all, or was not deployed at
% 0.6 m water depth)

UVAdates_good = UVAdates_all([2:10,13:15,20:37,41:45]);

% calculate daily integrated fluxes using same convention as NOAA ESRL, see
% http://esrl.noaa.gov/gmd/grad/antuv/docs/netOps/CHAPTER4.PDF, p. 4-25

% preallocate destination matrix

UVA_daily_dose_0_6m_subsurface_PAL1314_kJ_m2 = nan(length(min(UVAdates_all):max(UVAdates_all)),2);
UVA_daily_dose_0_6m_subsurface_PAL1314_kJ_m2(:,1) = datenum(min(UVAdates_all):max(UVAdates_all));

% make calculations with center of each integration period at local noon
% for Palmer, NOAA protocol (link above) defines this as approx. 1600 UTC
%
% note that the JAZ data J.R.C. recorded during the 2013-2014 field season
% at Palmer have timestamps in UTC, *not* local time

% identities

W_per_uW = 1/1000000;
J_per_kJ = 1/1000;
cm2_per_m2 = 10000;

for i=1:size(UVA_daily_dose_0_6m_subsurface_PAL1314_kJ_m2,1)
    if ismember(UVA_daily_dose_0_6m_subsurface_PAL1314_kJ_m2(i,1),datenum(UVAdates_good))
        % define the local noon for this date in MATLAB Julian format
        thisLocalNoon = datenum(UVA_daily_dose_0_6m_subsurface_PAL1314_kJ_m2(i,1)+datenum(0,0,0,16,0,0));
        % select the data within the local noon +/- 12 hr window
        theseUVAData = ...
            UVA_flux_PAL1314_uW_cm2(...
            UVA_flux_PAL1314_uW_cm2(:,1)>(thisLocalNoon-datenum(0,0,0,12,0,0)) & ...
            UVA_flux_PAL1314_uW_cm2(:,1)<(thisLocalNoon+datenum(0,0,0,12,0,0)),:);
        % calculate the time integrated-dose, in kJ/m2
        UVA_daily_dose_0_6m_subsurface_PAL1314_kJ_m2(i,2) = ...
            trapz(theseUVAData(:,1)*24*60*60,theseUVAData(:,2))*...
            W_per_uW*J_per_kJ*cm2_per_m2;
    end
end

% save .mat file

save('Daily_int_UVA_dose_0.6m_subsurface_PAL1314_kJ_m2.mat','UVA_daily_dose_0_6m_subsurface_PAL1314_kJ_m2')

% export UVB doses to .csv; have to use dlmwrite instead of csvwrite to get
% proper precision

dlmwrite('Daily_int_UVA_dose_0.6m_subsurface_PAL1314_kJ_m2.csv', UVA_daily_dose_0_6m_subsurface_PAL1314_kJ_m2, 'delimiter', ',', 'precision', 8); 

%% For time-series data on a given day, e.g., 14 Dec 13

subset_ind=find(UVB_flux_PAL1314_uW_cm2(:,1)>=datenum(2013,12,14)...
    & UVB_flux_PAL1314_uW_cm2(:,1)<datenum(2013,12,15));
UVB_14Dec=UVB_flux_PAL1314_uW_cm2(subset_ind,:);
times=UVB_14Dec(:,1);

subset_ind=find(Insitu_spectra_PAL1314_uW_cm2_good(:,1)>=datenum(2013,12,14)...
    & Insitu_spectra_PAL1314_uW_cm2_good(:,1)<datenum(2013,12,15));
specdata_14Dec=Insitu_spectra_PAL1314_uW_cm2_good(subset_ind,:);
times=specdata_14Dec(:,1);

specdata_14Dec_UVB=specdata_14Dec(:,ind_UVB+4);
specdata_14Dec_UVA=specdata_14Dec(:,ind_UVA+4);
% Offset because there are 4 columns of nonspectral data in the
% Insitu_spectra_PAL1314_uW_cm2_good matrix

% Generate integrated spectral dosages (in uW/cm2) for UVA, UVB at each JAZ timepoint

UVB_dose=nan(length(times),1);
UVA_dose=nan(length(times),1);

for i=1:length(times)
    UVB_dose(i,1) = trapz(UVB_wavelengths,specdata_14Dec_UVB(i,:));
    UVA_dose(i,1) = trapz(UVA_wavelengths,specdata_14Dec_UVA(i,:));
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
