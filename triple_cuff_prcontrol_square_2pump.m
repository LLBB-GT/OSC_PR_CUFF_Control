clear
clc

% Initialize the NIDAQ for data acquisition
ni_daq = daq.getDevices;
dev_ID = ni_daq.ID;

% Define the input channel for cuff 1
pr_sensor_1 = daq.createSession('ni'); 
addAnalogInputChannel(pr_sensor_1, dev_ID, 4, 'Voltage');
% Define the output channel for switching solenoid connected to cuff 1
valve_switch_1 = daq.createSession('ni');
addDigitalChannel(valve_switch_1, dev_ID, 'Port0/Line1','OutputOnly');

% Define the input channel for cuff 2
pr_sensor_2 = daq.createSession('ni');
addAnalogInputChannel(pr_sensor_2, dev_ID, 5, 'Voltage');
% Define the output channel for switching solenoid connected to cuff 2
valve_switch_2 = daq.createSession('ni');
addDigitalChannel(valve_switch_2, dev_ID, 'Port0/Line2','OutputOnly');

% Define the input channel for cuff 3
pr_sensor_3 = daq.createSession('ni');
addAnalogInputChannel(pr_sensor_3, dev_ID, 6, 'Voltage');
% Define the output channel for switching solenoid connected to cuff 3
valve_switch_3 = daq.createSession('ni');
addDigitalChannel(valve_switch_3, dev_ID, 'Port0/Line3','OutputOnly');

% Define the output channels for switching the solenoids
% connected to the syringe pumps
% For pump 1
valve_switch_4 = daq.createSession('ni');
addDigitalChannel(valve_switch_4, dev_ID, 'Port0/Line6','OutputOnly');
% For pump 2
valve_switch_5 = daq.createSession('ni');
addDigitalChannel(valve_switch_5, dev_ID, 'Port0/Line7','OutputOnly');

% Initialize syringe pump 1
s_1 = serial('COM5');
set(s_1,'BaudRate',115200);
fopen(s_1);
% Initialize syringe pump 2
s_2 = serial('COM6');
set(s_2,'BaudRate',115200);
fopen(s_2);

%% Move syringe to the starting position
start_pos = 11;   % Starting position, calculated in cm from the zero position
start_com = ['mov abs ' num2str(start_pos) ' cm 100%'];   % Command to move the pump to the starting position
fprintf(s_1,'?rep pos');   % Check the current position of pump 1
act_pos_1 = fscanf(s_1);
fprintf(s_2,'?rep pos');   % Check the current position of pump 2
act_pos_2 = fscanf(s_2);

% Get the current position of pump 1 in mm
act_pos_1(1:3) = [];
act_del_1 = findstr(act_pos_1,' mm');
act_pos_1(act_del_1:length(act_pos_1)) = [];
act_start_pos_1 = str2double(act_pos_1);   % Current position of pump 1 in mm
% Get the starting position of pump 2 in mm
act_pos_2(1:3) = [];
act_del_2 = findstr(act_pos_2,' mm');
act_pos_2(act_del_2:length(act_pos_2)) = [];
act_start_pos_2 = str2double(act_pos_2);   % Current position of pump 1 in mm

% Set solenoids to open to air
off_state = 0;

outputSingleScan(valve_switch_1, off_state);
outputSingleScan(valve_switch_2, off_state);
outputSingleScan(valve_switch_3, off_state);

% Start moving the pumps to the desired starting position
start = tic;
while abs(act_start_pos_1/10 - start_pos) > 0.05 ||...
        abs(act_start_pos_2/10 - start_pos) > 0.05   % Move the syringe if the desired starting position is 0.5 mm away from the actual position
    % Move syringe pump 1
    fprintf(s_1,start_com)
    pause(0.25)
    fprintf(s_1,'?rep pos');
    act_pos_1 = fscanf(s_1);
    act_pos_1(1:3) = [];
    act_del_1 = findstr(act_pos_1,' mm');
    act_pos_1(act_del_1:length(act_pos_1)) = [];
    act_start_pos_num = str2double(act_pos_1);
    if isnan(act_start_pos_num) == 0
        act_start_pos_1 = act_start_pos_num;   % Updated actual position for pump 1
    end

    % Move syringe pump 2
    fprintf(s_2,start_com)
    pause(0.25)
    fprintf(s_2,'?rep pos');
    act_pos_2 = fscanf(s_2);
    act_pos_2(1:3) = [];
    act_del_2 = findstr(act_pos_2,' mm');
    act_pos_2(act_del_2:length(act_pos_2)) = [];
    act_start_pos_num = str2double(act_pos_2);
    if isnan(act_start_pos_num) == 0
        act_start_pos_2 = act_start_pos_num;   % Updated actual position for pump 2
    end
    
    pause(0.5)
    
    elapsed = toc(start);
    % Stop the pumps after 25 s
    if elapsed > 25
        break
    end
end

fprintf(s_1,'stop')
fprintf(s_2,'stop')

%% Solenoid diagnostics - use this section to test the solenoids

on_state = 1;
off_state = 0;

outputSingleScan(valve_switch_1, on_state);
pause(1)
outputSingleScan(valve_switch_1, off_state);
pause(1)

outputSingleScan(valve_switch_2, on_state);
pause(1)
outputSingleScan(valve_switch_2, off_state);
pause(1)

outputSingleScan(valve_switch_3, on_state);
pause(1)
outputSingleScan(valve_switch_3, off_state);
pause(1)

outputSingleScan(valve_switch_4, on_state);
pause(1)
outputSingleScan(valve_switch_4, off_state);
pause(1)

outputSingleScan(valve_switch_5, on_state);
pause(1)
outputSingleScan(valve_switch_5, off_state);

%% Send the pressure information

% Define the pressure waveforms
n_wave = 3;   % Number of cuffs
amp = 60;   % Amplitude in mmHg
freq = 0.05;   % Frquency of each cuff in Hz
runtime = 600;   % Total duration of the pressure wave in seconds
dt = 0.1;   % Time resolution of the pressure waveform
T_frac = 3;   % Fraction of time period (to set delay between the pressure cuffs)
td = round(1/(T_frac*freq),1);   % Time delay between the pressure cuffs
t = 0:dt:runtime-dt+td*(n_wave-1);   % Total time (s)
duty = 50;                       % Duty cycle of the waveform in %
% Initialize the waveforms
wave_1 = zeros(1,length(t)); 
wave_1(1,1:round(runtime/dt)) = amp/2*(1+square(2*pi*freq*t(1,1:runtime/dt),duty));
wave_2 = zeros(1,length(t)); 
wave_2(1,td/dt+1:round((runtime+td)/dt)) = amp/2*(1+square(2*pi*freq*t(1,1:runtime/dt),duty));
wave_3 = zeros(1,length(t));
wave_3(1,2*td/dt+1:round((runtime+2*td)/dt)) = amp/2*(1+square(2*pi*freq*t(1,1:runtime/dt),duty));

% Separate the signal into 1 minute windows for switching between the
% syringe pumps
run_min = runtime/60;   % Number of windows
wind = 60;   % Window duration in seconds
% Initialize the pressure data array
pr_data_1 = [];
pr_data_2 = [];
pr_data_3 = [];

% Start controlling the pumps to achieve the desired pressure waveform
tot_start = tic;
for a = 1:run_min
    % Select which syringe pump to run, while the other resets, by setting
    % the valve states for the solenoids connected to the syringe pumps
    if mod(a,2) == 1
        s_pump = s_1;
        s_prime = s_2;
        outputSingleScan(valve_switch_4, on_state);
        outputSingleScan(valve_switch_5, off_state);
    elseif mod(a,2) == 0
        s_pump = s_2;
        s_prime = s_1;
        outputSingleScan(valve_switch_4, off_state);
        outputSingleScan(valve_switch_5, on_state);
    end
    start_pos = 11;
    start_com = ['mov abs ' num2str(start_pos) ' cm 100%'];
    fprintf(s_prime,start_com)
    
    % Send the move signal to the syringe pump
    Fs = 1;   % Sampling rate for the pressure sensor
    % Select the bounds for the syringe movement
    max_right = 1;
    max_left = 11;
    
    % Define the on and off states of the solenoids
    on_state = 1;
    off_state = 0;
    % Get the starting pressures (in Volts) in the three pressure sensors
    pr_data_1 = [pr_data_1 inputSingleScan(pr_sensor_1)];
    pr_data_2 = [pr_data_2 inputSingleScan(pr_sensor_2)];
    pr_data_3 = [pr_data_3 inputSingleScan(pr_sensor_3)];
    err_thresh = 1;   % in mmHg
    
    % Calculate the pressure and move the syringe pump within the 1 min
    % time windows
    for i = wind*(a-1)/dt+2:wind*a*Fs/dt
        start = tic;

        % Select the states of the solenoids connected to the pressure cuffs
        count_1 = 0;
        count_2 = 0;
        count_3 = 0;
        if wave_1(1,floor(i/Fs)+1) == amp
            outputSingleScan(valve_switch_1, on_state);
            count_1 = count_1 + 1;
        else
            outputSingleScan(valve_switch_1, off_state);
        end
        if wave_2(1,floor(i/Fs)+1) == amp
            outputSingleScan(valve_switch_2, on_state);
            count_2 = count_2 + 1;
        else
            outputSingleScan(valve_switch_2, off_state);
        end
        if wave_3(1,floor(i/Fs)+1) == amp
            outputSingleScan(valve_switch_3, on_state);
            count_3 = count_3 + 1;
        else
            outputSingleScan(valve_switch_3, off_state);
        end

        % Calculate the pressure in sensor 1 and also the error compared to
        % the desired pressure
        pr_data_1 = [pr_data_1 inputSingleScan(pr_sensor_1)];
        pr_data_gauge_1 = 51.71484*(7.5*(pr_data_1-0.5)-15);   % Voltage converted to psi and then to mmHg
        err_1 = pr_data_gauge_1(1,i) - wave_1(1,floor(i/Fs)+1);
        err_abs_1 = abs(err_1);

        % Calculate the pressure in sensor 2 and also the error compared to
        % the desired pressure
        pr_data_2 = [pr_data_2 inputSingleScan(pr_sensor_2)];
        pr_data_gauge_2 = 51.71484*(7.5*(pr_data_2-0.5)-15);   % Voltage converted to psi and then to mmHg
        err_2 = pr_data_gauge_2(1,i) - wave_2(1,floor(i/Fs)+1);
        err_abs_2 = abs(err_2);

        % Calculate the pressure in sensor 3 and also the error compared to
        % the desired pressure
        pr_data_3 = [pr_data_3 inputSingleScan(pr_sensor_3)];
        pr_data_gauge_3 = 51.71484*(7.5*(pr_data_3-0.5)-15);   % Voltage converted to psi and then to mmHg
        err_3 = pr_data_gauge_3(1,i) - wave_3(1,floor(i/Fs)+1);
        err_abs_3 = abs(err_3);

        % Control the pressure cuffs by moving the syringe pump as a
        % function of the error in the pressure
        % Control pressure cuff 1
        if count_1 == 1
            if err_1 < -1*err_thresh
                com = ['mov abs ' num2str(max_right) ' cm 100%'];
                fprintf(s_pump,com)
            elseif err_1 > err_thresh
                com = ['mov abs ' num2str(max_left) ' cm 100%'];
                fprintf(s_pump,com)
            elseif err_1 > -1*err_thresh && err_1 < 0
                com = ['mov abs ' num2str(max_right) ' cm ' num2str(floor(err_abs_1*100)) '%'];
                fprintf(s_pump,com)
             elseif err_1 > 0 && err_1 < err_thresh
                com = ['mov abs ' num2str(max_left) ' cm ' num2str(floor(err_abs_1*100)) '%'];
                fprintf(s_pump,com)
            end
        end

        % Control pressure cuff 2
        if count_2 == 1
            if err_2 < -1*err_thresh
                com = ['mov abs ' num2str(max_right) ' cm 100%'];
                fprintf(s_pump,com)
            elseif err_2 > err_thresh
                com = ['mov abs ' num2str(max_left) ' cm 100%'];
                fprintf(s_pump,com)
            elseif err_2 > -1*err_thresh && err_2 < 0
                com = ['mov abs ' num2str(max_right) ' cm ' num2str(floor(err_abs_2*100)) '%'];
                fprintf(s_pump,com)
             elseif err_2 > 0 && err_2 < err_thresh
                com = ['mov abs ' num2str(max_left) ' cm ' num2str(floor(err_abs_2*100)) '%'];
                fprintf(s_pump,com)
            end
        end

        % Control pressure cuff 3
        if count_3 == 1
            if err_3 < -1*err_thresh
                com = ['mov abs ' num2str(max_right) ' cm 100%'];
                fprintf(s_pump,com)
            elseif err_3 > err_thresh
                com = ['mov abs ' num2str(max_left) ' cm 100%'];
                fprintf(s_pump,com)
            elseif err_3 > -1*err_thresh && err_3 < 0
                com = ['mov abs ' num2str(max_right) ' cm ' num2str(floor(err_abs_3*100)) '%'];
                fprintf(s_pump,com)
             elseif err_3 > 0 && err_3 < err_thresh
                com = ['mov abs ' num2str(max_left) ' cm ' num2str(floor(err_abs_3*100)) '%'];
                fprintf(s_pump,com)
            end
        end
        elapsed = toc(start);
        pause(dt/Fs - elapsed)
    end
    fprintf(s_prime,'stop')
end
tot_elapsed = toc(tot_start);
serialbreak(s_pump)
fprintf(s_pump,'stop')

outputSingleScan(valve_switch_1, off_state);
outputSingleScan(valve_switch_2, off_state);
outputSingleScan(valve_switch_3, off_state);
outputSingleScan(valve_switch_4, off_state);
outputSingleScan(valve_switch_5, off_state);

% Release the devices
fclose(s_1);
fclose(s_2);
fclose(s_pump);
fclose(s_prime);
release(pr_sensor_1)
release(pr_sensor_2)
release(pr_sensor_3)

 %% Plot in the axis provided
fig_pathname = 'C:\Users\LLBB_NIR\Desktop\Anish\01.23.19\';
set(0,'defaultfigurecolor',[1 1 1]) % Set figure background to white
pr_time = linspace(0,tot_elapsed,length(pr_data_gauge_1));
plot(pr_time,pr_data_gauge_1,'k')
hold on
plot(t,wave_1,'r')
plot(pr_time,pr_data_gauge_2,'k--')
hold on
plot(t,wave_2,'r--')
plot(pr_time,pr_data_gauge_3,'k-.')
hold on
plot(t,wave_3,'r-.')
set(gca,'FontSize',12)
legend('Pr_{actual} Cuff 1','Pr_{desired} Cuff 1','Pr_{actual} Cuff 2',...
    'Pr_{desired} Cuff 2','Pr_{actual} Cuff 3','Pr_{desired} Cuff 3')
xlabel('Time (sec)','FontSize',18)
ylabel('Pressure (mmHg)','FontSize',18)
hold off
fig_filename = 'Rat_10ul_f02_1sec_40mmHg_Test_2.png';
fig_savename = [fig_pathname fig_filename];
% export_fig(fig_savename, '-m5')