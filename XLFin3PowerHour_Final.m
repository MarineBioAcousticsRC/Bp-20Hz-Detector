% Automatic acoustic energy detector for fin whale 20 Hz calls
%
% Edited by: Sofia Aniceto, Postdoctoral fellow at the Norwegian University
% of Science and Technology (NTNU)
% Project manager: Ana Sirovic
% 
% Calculates the difference in acoustic energy between the signal (22 Hz)
% and background noise levels (average energy between 10 and 34 Hz) from a 
% long-term spectral average (LTSA) in Triton.
% 
% Manual or automatic detection of fin whale 20 Hz calls is difficult when
% many animals are calling in an area at the same time. The calls overlap 
% and it can be impossible to determine where one call ends and another 
% begins. Because the calls appear as a band of energy, we detect average 
% fin whale acoustic energy instead of individual calls.
%
% Process: https://urldefense.com/v3/__https://github.com/MarineBioAcousticsRC/Bp-20Hz-Detector__;!!Mih3wA!HaJjuAmXPiaBJIA7Kz0wwhYNXSQejsvK22fE9kI2hLQ50g0PNifzywjDhMxFx4-OKwJ5LWHz7EpdVPFOCrjzHg$ 
% 1 - Open triton through matlab
% 2 - Open ltsa file
% 3 - Change datetime of the first ltsa loaded (if calculating this for several) to make sure it's at 00:00:00 HMS
% 4 - Change plot length to 1h
% 5 - Run this script
%
% IMPORTANT!!! DON'T FORGET TO CHANGE THE LTSA PLOT LENGTH IN TRITON!!!
% OTHERWISE IT PRODUCES GAPS IN THE OUTPUT
%
% Output: rather than the XML mentioned in the wiki above, here it's csv
% The script produces a plot for hourly estimates, with raw data in the csv
%

clc
clearvars

%% Load the data
global PARAMS                                            
tic

%% Directories
% Data directory 
xl_dir = '\\home.ansatt.ntnu.no\anasan\Documents\SOCAL data analyses\Bp_20Hz\SOCAL_E_71_preamp966\';
new_filename = PARAMS.ltsa.infile;
xl_filename = new_filename(1:strfind(new_filename,'.ltsa')-1); % remove ltsa from filename
xl_out= strcat(xl_dir,xl_filename,'20finDPI_hour.csv'); % Output csv file

%% Read LTSA
% Opening LTSA file to be able to continue reading from it
fid = fopen([PARAMS.ltsa.inpath,PARAMS.ltsa.infile],'r');

%% Get information and wrangle ltsa data
LTSAres_time = num2str(PARAMS.ltsa.tave); % LTSA time resolution (s) 
LTSAres_freq = num2str(PARAMS.ltsa.dfreq); % LTSA frequency resolution (Hz)
LTSA_loaded_sec = PARAMS.ltsa.tseg.sec; % Total duration of the loaded LTSA (s)
LTSA_loaded_dur = PARAMS.ltsa.tseg.hr; % Total loaded duration of the loaded LTSA (h)
LTSA_duration_timestamps = PARAMS.ltsa.dur; % Total duration with timestamps of the LTSA
LTSA_duration_time = (33310/3600*75)/24; % duration of entire ltsa in days

%LTSA_all_bin = size(LTSA_duration_timestamps,2)*15; %(75/5) %how many 5s bins are in the 1h ltsa 
%3330*15; %nr of 5s bins in the 33310 of 75s

% Each 2-h file (1440 samples) has only an absolute time vector (from 0
% to 2h) saved
% The HARP writes data to the disk every 75s, and PARAMS.ltsa.dnumStart
% is written every 75s. So we need the 5s power estimates to match the
% timestamps split in 75s
HARP_samples = 75;

% Begin time of the LTSA
ltsaplotstart = PARAMS.ltsa.plotStartRawIndex; % index
firsttime = PARAMS.ltsa.dnumStart(ltsaplotstart); % dnum means datenum in days (from sm_read_ltsahead.m) 
disp(['firsttime: ', dbSerialDateToISO8601(firsttime)]) 

% Time split
% PARAMS.ltsa.tseg.hr is initial window time segment duration in hours
% by default, it opens 2h but let's load 1h (3600 seconds) of LTSA
% when the LTSA is made (set in mk_ltsa.m)
% PARAMS.ltsa.tave is the averaging time [5 seconds] for the LTSA (set in init_ltsaparams.m)
nbin = floor(((LTSA_loaded_dur) * 60 *60 ) / PARAMS.ltsa.tave); % number of 5-s time bins in one hour ltsa

% Skip the hearder and start of ltsa because of self noise

% time bin number at start of plot within rawfile (index)
skip = PARAMS.ltsa.byteloc(PARAMS.ltsa.plotStartRawIndex) + ...
   (PARAMS.ltsa.plotStartBin - 1) * PARAMS.ltsa.nf + ...
    PARAMS.ltsa.nf * nbin;
fseek(fid,skip,-1);    % skip over header + other data read in first LTSA


%% Calibration
% Load appropriate transfer function to apply to data
tffilename =  '\\home.ansatt.ntnu.no\\anasan\\Documents\\SOCAL data analyses\\Bp_20Hz\\SOCAL_E_71_preamp966\\transfer functions\\966_210115_A_HARP.tf';
fidtf1 = fopen(tffilename,'r');
[A1,~] = fscanf(fidtf1,'%f %f',[2,inf]);
freqvec = 1:1001;
tf_HARP = interp1(A1(1,1:60),A1(2,1:60),freqvec,'linear','extrap');


% Apply the transfer function to the ltsa duration 
pwr = PARAMS.ltsa.pwr; %PARAMS.lta.pwr is a matrix [1001,720], 720 is the nr of 5 second bins in one hour of data (3600/5)
%to reformat this to represent only one hour, we need to cut the pwr vector
%to reflect the first 1h of the ltsa but keep all the original frequencies
pwr = pwr + tf_HARP';
clear tffilename fidtf1 A1 count1 freqvec



%% Detector parameters

callfreq = 22;
nfreq1 = 10;
nfreq2 = 34;

% Record start of effort
effort(1) = firsttime+datenum([2000 0 0 0 0 0]);
effStart = dbSerialDateToISO8601(effort(1));
tc = find(PARAMS.ltsa.dnumStart==firsttime); %PARAMS.ltsa.dnumStart has the timestamps of every 75s bins. tc is an index for where the start of the ltsa with 75s timestamps meets the start of the loaded ltsa


%% Run
% This next while loop calculates the power for the 1h chunks that are loaded in
% the ltsa (so, instead of using the 2h that are loaded by default)
% PARAMS.ltsa.dnumStart is the starting time [dnum => datenum in days] for each ltsa raw

% Variable initialization for the loop
% create empty variables
ttime = []; totalpwer = []; pwr_avg= []; indx = []; newvec = []; scoreVal = []; hourVal = []; StartTime = [];


% This goes through each file in fid - each file lasts 2h
while ~feof(fid)
    
    %% Power
    % Calculate the SNR between freq in fin band and adjacent noise
    % (averaged from freq band above and below calls and assumed linear)
    pwr_avg = pwr(callfreq,:)-((pwr(nfreq1,:)+ pwr(nfreq2,:))/2); % [1x720]
    % find all times when the SNR is negative and fix to 0 
    pwr_avg(pwr_avg<0)=0;
    
    % because each file is transfered to the disk every 75s, we need to
    % estimate the nr of 5s bins in those 75s. So, we divide the total nr of samples by 
    % the nr of 5s bins in 75s resolution (15)
    stepss =  size(pwr,2)/(HARP_samples/PARAMS.ltsa.tave); % This how many samples are used for averaging - because the timestamps are stored in 75s bins (HARP_samples) and we have power estimates every 5s

    % Now we need to create a matrix that stores the information
    %so we resize the pwr_avg [1,720] to represent the stepss in
    %calculation (48)
    %by estimating the nr of 5s bins in the 75s samples (75/5=15)
    newvec = reshape(pwr_avg,[15,size(pwr_avg,2)/15]); % Reshapes the data from 1 x 720 into a 15 x 48 matrix, 
    newvec = mean(newvec,1); %then takes the average of each line do total power is 1 x 48 
  
    totalpwer = [totalpwer newvec]; % Save the vector in totalpwer [1 x 48]
    
    %% Time
    % assign time stamps to each 75s chunck
    % tc is where the start of the loaded ltsa meets the timestamp marked in the 75s
    % stepss is vector for the steps needed to consider the 5s averages
    % PARAMS.ltsa.dnumStart is starting time [dnum => datenum in days] for each ltsa raw
    if (tc+stepss)<=size(PARAMS.ltsa.dnumStart,2) %if within the entire duration of the ltsa (with 75s timestamps)
        ttime = [ttime; PARAMS.ltsa.dnumStart(tc:tc+stepss-1)'];
    else
        ttime = [ttime; PARAMS.ltsa.dnumStart(tc:end)'];
    end
    tc = tc+stepss;
 
    ttimevec = datevec(ttime); % [YY, MM, DD, HH, MM, SS] one day of timestamps %should be [48,6]

    % When it goes into a new hour within the first day, we average it all out
    if ttimevec(1,4) ~= ttimevec(size(ttimevec,1),4) % checking if we're in a new hour
        % Average the power and add datetime stamp
        indx = find(ttimevec(:,4)==ttimevec(1,4)); %should be [48,1]
        hourlypwr = mean(totalpwer(indx)); %should have one value
        scoreVal = [scoreVal hourlypwr]; %should be one value
        timestamp = datenum(ttimevec(1,1),ttimevec(1,2),ttimevec(1,3),ttimevec(1,4),0,0); %should be one value
        
        hourVal = [hourVal timestamp+datenum([2000 0 0 0 0 0])]; %should be one value
        
        % need to save data that are not in that day/hour
        newpwr = totalpwer(indx(end)+1:end);
        newtime = ttime(indx(end)+1:end);
        
        %empty them before re-writing
        ttime = []; totalpwer = []; pwr_avg= []; indx = []; pwr = []; newvec = [];
       
        totalpwer = newpwr;
        ttime = newtime;
        
        %empty new time and newpwr
        newtime = []; newpwr = [];
        
        %then write the detection
        %includes dailypwr; timestamp+datenum([2000 0 0 0 0 0]);
        dtime = timestamp+datenum([2000 0 0 0 0 0]);
        %check that the beginning of effort falls before the timestamp
        startISO = dbSerialDateToISO8601(dtime);
        if dtime<effort(1)
            startISO = effStart;
        end
        score = hourlypwr;
    end
    
    % Read next chunk data
    pwr = fread(fid,[PARAMS.ltsa.nf,nbin],'int8');   
    
    % Apply the transfer function
    pwr = pwr + tf_HARP';
    
end


hourlypwr = mean(totalpwer);
scoreVal = [scoreVal hourlypwr];
ttimevec = datevec(ttime);
timestamp = datenum(ttimevec(1,1),ttimevec(1,2),ttimevec(1,3));

hourVal = [hourVal timestamp+datenum([2000 0 0 0 0 0])];
dtime = timestamp+datenum([2000 0 0 0 0 0]);
startISO = dbSerialDateToISO8601(dtime);
score = hourlypwr;
        
%record end of effort
effort(2) = PARAMS.ltsa.dnumStart(end);
effort(2) = effort(2)+datenum([2000 0 0 0 0 0]);
effEnd = dbSerialDateToISO8601(effort(2)); 

datetime = datetime(datestr(hourVal));

% datevec(effort) %
t = toc;
disp(' ')
disp(['FinPowerDetect Execution time = ',num2str(t),' seconds'])

% Plot
Times = cell2mat(dbSerialDateToISO8601(hourVal'));
disp(['Begin hour: ', Times(1,:),' // End hour: ', Times(end,:)])

figure;
%plot(hourVal,scoreVal);
plot(datetime,scoreVal);
%plot(hourVal,scoreVal,'*');
xlabel('Day');
ylabel('Fin acoustic power index');

% csv write
%csvwrite(xl_out,[hourVal' scoreVal'])
tst_data = table(datetime, hourVal', scoreVal'); 
writetable(tst_data, xl_out) 




