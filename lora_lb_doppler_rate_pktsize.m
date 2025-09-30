%% SATELLITE LINK BUDGET ANALYSIS
% Juan Del Pino Mena
% June 2025
%
% Link budget estimator for LEO satellites.
% This script computes the Doppler Rate for several spreading factors and orbit altitudes,
% sweeping the size of the LoRa packet.
%
% This program requires MatLab >= R2023a and the satellite communications toolbox.


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Script config

close all;
clearvars;

show3Dviewer = false;
printOptional = true;

warning('off', 'MATLAB:hgexport:NotSupportedInFutureRelease');
warning('off', 'satcom:p618PropagationLosses:InvalidXPDFrequency');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants

k_boltzmann = 1.380649e-23;  % [J/K] Boltzmann's constant
T_0 = 290; % [K] room temperature for noise calculus, usually 290 K
R_E = 6371;  % [km] Earth's average radius


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Configuration parameters

% ----------------------------------------------------------------------------------------
% Common transmission parameters

freq_Hz = 870e6;  % [Hz] Transmission frequency
bandwidth_Hz = 125e3;  % [Hz] Transmission bandwidth
symbolperiod_s = 0.03277;  % [s] Symbol period. USE 0.03277 AS MINIMUM, NOT LOWER!!

symbolrate_Sps = 1/symbolperiod_s;  % [Sps] Symbol rate (LoRa: Ts = 2^SF / BW)
bitrate_bps = 293;  % [b/s] Transmission bitrate

% ----------------------------------------------------------------------------------------
% Transceiver limits (only used for plots)

rx_limit_sensitivity_gs = -137;  % [dBm] The minimum admitted signal power in the gs
rx_limit_cnr_gs = -20;  % [dB] The minimum SNR / CNR the gs transceiver admits

rx_limit_sensitivity_sat = -140;  % [dBm] The minimum admitted signal power in the sat
rx_limit_cnr_sat = -20;  % [dB] The minimum SNR / CNR the sat transceiver admits

% ----------------------------------------------------------------------------------------
% Ground Station

gs_name = "UPV GS";
gs_lat_N = 39.47943;  % North latitude
gs_lon_E = -0.34230;  % East longitude
gs_altitude = 50;  % [m] Altitude of the GS
elevation_min_angle = 5;  % [degrees] Minimum angle of elevation in GS (defines access intervals)

% Power, gain and losses

tx_power_gs_dBm = 14;  % [dBm] Transmission power in the Ground Station tx
tx_ant_gain_gs_dB = 2;  % [dBi] Antenna Gain in the GS
tx_loss_gs_dB = 1;  % [dB] TX system power loss in the ground station
rx_loss_gs_dB = 3;  % [dB] RX system power loss in the ground station

% Noise

rx_nf_gs_dB = 6;  % [dB] System noise figure in RX in the ground station

ant_ambient_temp_gs_K = 300;  % [K] Antenna ambient temperature in the ground station
ant_noise_temp_gs_K = 300;  % [K] Antenna noise temperature in the ground station
% The equivalent noise temperature seen at the receive output of the antenna. Describes
% how much noise an antenna produces in a given environment.

% G/T (Gain-to-Noise-Temperature Ratio, GNTR) relates the receive antenna gain and the
% system noise temperature.
% G/T = GR − Nf − 10 * log10( T0 + (Ta − T0) * 10 ^ (−0.1 * Nf) )
% GR is the receive antenna gain in dBi
% Nf is the noise figure in dB
% T0 is the ambient temperature in degrees Kelvin
% Ta is the antenna temperature in degrees Kelvin
% https://www.mathworks.com/help/satcom/ug/nb-iot-ntn-link-budget-analysis.html

rx_GNTR_gs_dB_K = tx_ant_gain_gs_dB - rx_nf_gs_dB - 10 * log10(ant_ambient_temp_gs_K + ...
    (ant_noise_temp_gs_K - ant_ambient_temp_gs_K) * 10^(-0.1 * rx_nf_gs_dB));  % [dB/K]

% ----------------------------------------------------------------------------------------
% Satellite

% Power, gain and losses

tx_power_sat_dBm = 33;  % [dBm] Transmission power in the satellite tx
tx_ant_gain_sat_dB = 2;  % [dBi] Antenna Gain in the sat
tx_loss_sat_dB = 2;  % [dB] System power loss in TX in the satellite
rx_loss_sat_dB = 1.5;  % [dB] System power loss in RX in the satellite

% Noise

rx_nf_sat_dB = 2;  % [dB] System noise figure in RX in the satellite

ant_ambient_temp_sat_K = 350;  % [K] Antenna ambient temperature in the satellite
ant_noise_temp_sat_K = 350;  % [K] Antenna noise temperature in the satellite

% [dB/K] G/T in the sat rx
rx_GNTR_sat_dB_K = tx_ant_gain_sat_dB - rx_nf_sat_dB ...
    - 10 * log10(ant_ambient_temp_sat_K ...
    + (ant_noise_temp_sat_K - ant_ambient_temp_sat_K) * 10^(-0.1 * rx_nf_sat_dB));

% ----------------------------------------------------------------------------------------
% Simulation time, stop time and step

dl = datetime('25-Sep-2024 02:25:00', 'Format', 'yyyy-MM-dd HH:mm:ss', TimeZone="UTC");
dr = dl + minutes(30); % datetime('25-Sep-2024 00:20:00', 'Format', 'yyyy-MM-dd HH:mm:ss', TimeZone="UTC");
sampleTime = 0.1;  % [s] Simulation step (a small step is necessary for doppler rate calc)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create scenario

scenario = satelliteScenario(dl, dr, sampleTime);  % Satellite scenario


% ----------------------------------------------------------------------------------------
% Satellite
altitude = 300;  % km
semiMajorAxis = (R_E + altitude) * 1e3; % [m]
eccentricity = 0;
inclination = 50;  % [deg]
rightAscensionOfAscendingNode = 0;  % [deg]
argumentOfPeriapsis = 0;  % [deg]
trueAnomaly = 0;  % [deg]

sat = satellite(scenario, ...
                semiMajorAxis, ...
                eccentricity, ...
                inclination, ...
                rightAscensionOfAscendingNode,...
                argumentOfPeriapsis, ...
                trueAnomaly, ...
                Name="sat");


% ------------------------------------------------------------------------------
% Ground Station
gs = groundStation(scenario, Name=gs_name, Latitude=gs_lat_N, Longitude=gs_lon_E, ...
                   Altitude=gs_altitude, MinElevationAngle=elevation_min_angle);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add TX and RX to the satellite and the Ground Station

% ------------------------------------------------------------------------------
% Ground Station transceivers

gs_antenna = arrayConfig(Size=[1, 1]);  % Isotropic antenna from phasedArray T.Box

% Transmitter
gs_transmitter = transmitter( ...
    gs, ...  % parent element
    Antenna = gs_antenna, ...  % Antenna element
    Frequency = freq_Hz, ...  % Center frequency (Hz)
    Power = tx_power_gs_dBm - 30, ...  % TX power (dBW)
    BitRate = bitrate_bps * 1e-6, ...  % Bitrate (Mbps)
    SystemLoss = tx_loss_gs_dB, ...  % Equipment losses (dB)
    MountingAngles = [0; 0; 0], ...  % Orientation (degrees)
    Name = "GS TX" ...  % Name
);

% Receiver
gs_receiver = receiver( ...
    gs, ...  % parent element
    Antenna = gs_antenna, ...  % Antenna element
    MountingAngles = [0; 0; 0], ...  % Orientation (degrees)
    SystemLoss = rx_loss_gs_dB, ...  % Equipment losses (dB)
    Name = "GS RX" ...  % Name
);

% ------------------------------------------------------------------------------
% Satellite transceivers

sat_antenna = arrayConfig(Size=[1, 1]);  % Isotropic antenna

% Transmitter
sat_transmitter = transmitter( ...
    sat, ...  % parent element
    Antenna = sat_antenna, ...  % Antenna element
    Frequency = freq_Hz, ...  % Center frequency (Hz)
    Power = tx_power_sat_dBm - 30, ...  % TX power (dBW)
    BitRate = bitrate_bps * 1e-6, ...  % Bitrate (Mbps)
    SystemLoss = tx_loss_sat_dB, ...  % Equipment losses (dB)
    MountingAngles = [0; 0; 0], ...  % Orientation (degrees)
    Name = "SAT TX" ...  % Name
);

% Receiver
sat_receiver = receiver( ...
    sat, ...  % parent element
    Antenna = sat_antenna, ...  % Antenna element
    MountingAngles = [0; 0; 0], ...  % Orientation (degrees)
    SystemLoss = rx_loss_sat_dB, ...  % Equipment losses (dB)
    Name = "SAT RX" ...  % Name
);




%% Configure loop

loop_iter = 0;

payload_n_bytes_loop = 32:64:256;

npreamble = 8; % Example value for npreamble
DE = 1;         % Example value for DE (Low Data Rate Optimization, 0 or 1)
CR = 1;         % Example value for CR (Coding Rate, 1 to 4)
BW = 125e3;     % Example value for BW (Bandwidth in Hz)
H = 0;

SF = [12, 11, 10];


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEGIN FOR LOOP
for payload_n_bytes=payload_n_bytes_loop
    loop_iter = loop_iter + 1;

    PL = payload_n_bytes;

    sym12_t = 0.03277;  % [s] how long a symbol of sf12 lasts
    sym11_t = sym12_t/2;  % [s] how long a symbol of sf11 lasts
    sym10_t = sym12_t/4;  % [s] how long a symbol of sf10 lasts

    Tsym = (2.^SF) ./ BW;
    Tpreamble = (npreamble + 4.25) .* Tsym;
    payloadSymbNb = 8 + max(ceil((8 .* PL - 4 .* SF + 28 + 16 - 20 .* H) ./ (4 .* (SF - 2 .* DE))) .* (CR + 4), 0);
    Tpayload = payloadSymbNb .* Tsym;
    Tpacket = Tpreamble + Tpayload;

    Npacket = ceil(Tpacket./sampleTime);  % how many samples does the packet last

    pkt12_t = Tpacket(1);  % how many samples does the packet last
    pkt11_t = Tpacket(2);
    pkt10_t = Tpacket(3);

    pkt12_n = Npacket(1);  % how many samples does the packet last
    pkt11_n = Npacket(2);
    pkt10_n = Npacket(3);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Latency and Doppler Shift
    % These functions are only available on MatLab R2023a onwards
    % Dynamic - computes latency and doppler shift over the complete simulation time

    [doppler_fshift, doppler_time, doppler_info] = dopplershift(sat, gs, Frequency=freq_Hz);

    doppler_fshift_loop(loop_iter, :) = doppler_fshift;
    doppler_time_loop(loop_iter, :) = doppler_time;
    doppler_info_loop(loop_iter, :) = doppler_info;

    % -----------------------
    % doppler rate, old receivers (SX12xx)

    % https://www.semtech.com/design-support/lora-calculator
    % 255 bytes, 125 khz, 868 MHz

    pkt12_ddiff = zeros(1, length(doppler_time) - pkt12_n);  % alloc
    pkt11_ddiff = zeros(1, length(doppler_time) - pkt11_n);
    pkt10_ddiff = zeros(1, length(doppler_time) - pkt10_n);

    % calculate the doppler rate (substract samples pkt12_n samples apart)

    for i = 1:(length(doppler_time) - pkt12_n)
        pkt12_ddiff(i) = doppler_fshift(i+pkt12_n) - doppler_fshift(i);
    end
    pkt12_ddiff = abs(pkt12_ddiff);

    for i = 1:(length(doppler_time) - pkt11_n)
        pkt11_ddiff(i) = doppler_fshift(i+pkt11_n) - doppler_fshift(i);
    end
    pkt11_ddiff = abs(pkt11_ddiff);

    for i = 1:(length(doppler_time) - pkt10_n)
        pkt10_ddiff(i) = doppler_fshift(i+pkt10_n) - doppler_fshift(i);
    end
    pkt10_ddiff = abs(pkt10_ddiff);


    figure(1);

    subplot(3, 2, 1);
    plot(doppler_time(1:end-pkt12_n), pkt12_ddiff);
    hold on;

    subplot(3, 2, 3);
    plot(doppler_time(1:end-pkt11_n), pkt11_ddiff);
    hold on;

    subplot(3, 2, 5);
    plot(doppler_time(1:end-pkt10_n), pkt10_ddiff);
    hold on;

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END FOR LOOP
%% GRAPHS


dl_g = dl + minutes(11);  % limits
dr_g = dr - minutes(12);

orbits_labels_sf12 = {};
orbits_labels_sf11 = {};
orbits_labels_sf10 = {};

i = 0;
for payload_n_bytes=payload_n_bytes_loop
    i = i + 1;

    SF=12;
    PL=payload_n_bytes;
    Tsym = (2.^SF) ./ BW;
    Tpreamble = (npreamble + 4.25) .* Tsym;
    payloadSymbNb = 8 + max(ceil((8 .* PL - 4 .* SF + 28 + 16 - 20 .* H) ./ (4 .* (SF - 2 .* DE))) .* (CR + 4), 0);
    Tpayload = payloadSymbNb .* Tsym;
    Tpacket = Tpreamble + Tpayload;

    orbits_labels_sf12{i} = sprintf("%d B, %1.2f s", payload_n_bytes, Tpacket);
end
i = 0;
for payload_n_bytes=payload_n_bytes_loop
    i = i + 1;

    SF=11;
    PL=payload_n_bytes;
    Tsym = (2.^SF) ./ BW;
    Tpreamble = (npreamble + 4.25) .* Tsym;
    payloadSymbNb = 8 + max(ceil((8 .* PL - 4 .* SF + 28 + 16 - 20 .* H) ./ (4 .* (SF - 2 .* DE))) .* (CR + 4), 0);
    Tpayload = payloadSymbNb .* Tsym;
    Tpacket = Tpreamble + Tpayload;

    orbits_labels_sf11{i} = sprintf("%d B, %1.2f s", payload_n_bytes, Tpacket);
end
i = 0;
for payload_n_bytes=payload_n_bytes_loop
    i = i + 1;

    SF=10;
    PL=payload_n_bytes;
    Tsym = (2.^SF) ./ BW;
    Tpreamble = (npreamble + 4.25) .* Tsym;
    payloadSymbNb = 8 + max(ceil((8 .* PL - 4 .* SF + 28 + 16 - 20 .* H) ./ (4 .* (SF - 2 .* DE))) .* (CR + 4), 0);
    Tpayload = payloadSymbNb .* Tsym;
    Tpacket = Tpreamble + Tpayload;

    orbits_labels_sf10{i} = sprintf("%d B, %1.2f s", payload_n_bytes, Tpacket);
end

% ------------------------------------------------------------------------------
% doppler rate limits for LoRa SX12xx and SX1301, LDRO=on

pkt12_lim = floor((16 * 125e3) / (3 * 2^12));
pkt11_lim = floor((16 * 125e3) / (3 * 2^11));
pkt10_lim = floor((16 * 125e3) / (3 * 2^10));

pkt_ylim = 5000;
pkt_tick = 1000;

figure(1);

subplot(3, 2, 1);
grid on; grid minor;
title(sprintf("Doppler rate, SX12xx, SF12, h = %d km", altitude));
xlabel("Sim time");
xlim([dl_g dr_g]);
ylabel("Doppler rate (Hz/T_{pkt})");
ylim([0, pkt_ylim]);
yticks(0:pkt_tick:pkt_ylim);
yline(pkt12_lim, Color="#5e5e5e", LineStyle="--", Label=sprintf("%d Hz", pkt12_lim));

subplot(3, 2, 3);
grid on; grid minor;
title(sprintf("Doppler rate, SX12xx, SF11, h = %d km", altitude));
xlabel("Sim time");
xlim([dl_g dr_g]);
ylabel("Doppler rate (Hz/T_{pkt})");
ylim([0, pkt_ylim]);
yticks(0:pkt_tick:pkt_ylim);
yline(pkt11_lim, Color="#5e5e5e", LineStyle="--", Label=sprintf("%d Hz", pkt11_lim));

subplot(3, 2, 5);
grid on; grid minor;
title(sprintf("Doppler rate, SX12xx, SF10, h = %d km", altitude));
xlabel("Sim time");
xlim([dl_g dr_g]);
ylabel("Doppler rate (Hz/T_{pkt})");
ylim([0, pkt_ylim]);
yticks(0:pkt_tick:pkt_ylim);
yline(pkt10_lim, Color="#5e5e5e", LineStyle="--", Label=sprintf("%d Hz", pkt10_lim));


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create scenario

scenario = satelliteScenario(dl, dr, sampleTime);  % Satellite scenario


% ----------------------------------------------------------------------------------------
% Satellite
altitude = 450;  % km
semiMajorAxis = (R_E + altitude) * 1e3; % [m]
eccentricity = 0;
inclination = 50;  % [deg]
rightAscensionOfAscendingNode = 0;  % [deg]
argumentOfPeriapsis = 0;  % [deg]
trueAnomaly = 0;  % [deg]

sat = satellite(scenario, ...
                semiMajorAxis, ...
                eccentricity, ...
                inclination, ...
                rightAscensionOfAscendingNode,...
                argumentOfPeriapsis, ...
                trueAnomaly, ...
                Name="sat");


% ------------------------------------------------------------------------------
% Ground Station
gs = groundStation(scenario, Name=gs_name, Latitude=gs_lat_N, Longitude=gs_lon_E, ...
                   Altitude=gs_altitude, MinElevationAngle=elevation_min_angle);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add TX and RX to the satellite and the Ground Station

% ------------------------------------------------------------------------------
% Ground Station transceivers

gs_antenna = arrayConfig(Size=[1, 1]);  % Isotropic antenna from phasedArray T.Box

% Transmitter
gs_transmitter = transmitter( ...
    gs, ...  % parent element
    Antenna = gs_antenna, ...  % Antenna element
    Frequency = freq_Hz, ...  % Center frequency (Hz)
    Power = tx_power_gs_dBm - 30, ...  % TX power (dBW)
    BitRate = bitrate_bps * 1e-6, ...  % Bitrate (Mbps)
    SystemLoss = tx_loss_gs_dB, ...  % Equipment losses (dB)
    MountingAngles = [0; 0; 0], ...  % Orientation (degrees)
    Name = "GS TX" ...  % Name
);

% Receiver
gs_receiver = receiver( ...
    gs, ...  % parent element
    Antenna = gs_antenna, ...  % Antenna element
    MountingAngles = [0; 0; 0], ...  % Orientation (degrees)
    SystemLoss = rx_loss_gs_dB, ...  % Equipment losses (dB)
    Name = "GS RX" ...  % Name
);

% ------------------------------------------------------------------------------
% Satellite transceivers

sat_antenna = arrayConfig(Size=[1, 1]);  % Isotropic antenna

% Transmitter
sat_transmitter = transmitter( ...
    sat, ...  % parent element
    Antenna = sat_antenna, ...  % Antenna element
    Frequency = freq_Hz, ...  % Center frequency (Hz)
    Power = tx_power_sat_dBm - 30, ...  % TX power (dBW)
    BitRate = bitrate_bps * 1e-6, ...  % Bitrate (Mbps)
    SystemLoss = tx_loss_sat_dB, ...  % Equipment losses (dB)
    MountingAngles = [0; 0; 0], ...  % Orientation (degrees)
    Name = "SAT TX" ...  % Name
);

% Receiver
sat_receiver = receiver( ...
    sat, ...  % parent element
    Antenna = sat_antenna, ...  % Antenna element
    MountingAngles = [0; 0; 0], ...  % Orientation (degrees)
    SystemLoss = rx_loss_sat_dB, ...  % Equipment losses (dB)
    Name = "SAT RX" ...  % Name
);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Configure loop

loop_iter = 0;

payload_n_bytes_loop = 32:64:256;

npreamble = 8; % Example value for npreamble
DE = 1;         % Example value for DE (Low Data Rate Optimization, 0 or 1)
CR = 1;         % Example value for CR (Coding Rate, 1 to 4)
BW = 125e3;     % Example value for BW (Bandwidth in Hz)
H = 0;

SF = [12, 11, 10];


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEGIN FOR LOOP
for payload_n_bytes=payload_n_bytes_loop
    loop_iter = loop_iter + 1;

    PL = payload_n_bytes;

    sym12_t = 0.03277;  % [s] how long a symbol of sf12 lasts
    sym11_t = sym12_t/2;  % [s] how long a symbol of sf11 lasts
    sym10_t = sym12_t/4;  % [s] how long a symbol of sf10 lasts

    Tsym = (2.^SF) ./ BW;
    Tpreamble = (npreamble + 4.25) .* Tsym;
    payloadSymbNb = 8 + max(ceil((8 .* PL - 4 .* SF + 28 + 16 - 20 .* H) ./ (4 .* (SF - 2 .* DE))) .* (CR + 4), 0);
    Tpayload = payloadSymbNb .* Tsym;
    Tpacket = Tpreamble + Tpayload;

    Npacket = ceil(Tpacket./sampleTime);  % how many samples does the packet last

    pkt12_t = Tpacket(1);  % how many samples does the packet last
    pkt11_t = Tpacket(2);
    pkt10_t = Tpacket(3);

    pkt12_n = Npacket(1);  % how many samples does the packet last
    pkt11_n = Npacket(2);
    pkt10_n = Npacket(3);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Latency and Doppler Shift
    % These functions are only available on MatLab R2023a onwards
    % Dynamic - computes latency and doppler shift over the complete simulation time

    [doppler_fshift, doppler_time, doppler_info] = dopplershift(sat, gs, Frequency=freq_Hz);

    doppler_fshift_loop(loop_iter, :) = doppler_fshift;
    doppler_time_loop(loop_iter, :) = doppler_time;
    doppler_info_loop(loop_iter, :) = doppler_info;

    % -----------------------
    % doppler rate, old receivers (SX12xx)

    % https://www.semtech.com/design-support/lora-calculator
    % 255 bytes, 125 khz, 868 MHz

    pkt12_ddiff = zeros(1, length(doppler_time) - pkt12_n);  % alloc
    pkt11_ddiff = zeros(1, length(doppler_time) - pkt11_n);
    pkt10_ddiff = zeros(1, length(doppler_time) - pkt10_n);

    % calculate the doppler rate (substract samples pkt12_n samples apart)

    for i = 1:(length(doppler_time) - pkt12_n)
        pkt12_ddiff(i) = doppler_fshift(i+pkt12_n) - doppler_fshift(i);
    end
    pkt12_ddiff = abs(pkt12_ddiff);

    for i = 1:(length(doppler_time) - pkt11_n)
        pkt11_ddiff(i) = doppler_fshift(i+pkt11_n) - doppler_fshift(i);
    end
    pkt11_ddiff = abs(pkt11_ddiff);

    for i = 1:(length(doppler_time) - pkt10_n)
        pkt10_ddiff(i) = doppler_fshift(i+pkt10_n) - doppler_fshift(i);
    end
    pkt10_ddiff = abs(pkt10_ddiff);


    figure(1);

    subplot(3, 2, 2);
    plot(doppler_time(1:end-pkt12_n), pkt12_ddiff);
    hold on;

    subplot(3, 2, 4);
    plot(doppler_time(1:end-pkt11_n), pkt11_ddiff);
    hold on;

    subplot(3, 2, 6);
    plot(doppler_time(1:end-pkt10_n), pkt10_ddiff);
    hold on;

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END FOR LOOP
%% GRAPHS


% ------------------------------------------------------------------------------
% doppler rate limits for LoRa SX12xx and SX1301, LDRO=on

figure(1);


dl_g = dl + minutes(12);  % limits
dr_g = dr - minutes(10);

subplot(3, 2, 2);
grid on; grid minor;
title(sprintf("Doppler rate, SX12xx, SF12, h = %d km", altitude));
ylabel("Doppler rate (Hz/T_{pkt})");
xlabel("Sim time");
xlim([dl_g dr_g]);
ylim([0, pkt_ylim]);
yticks(0:pkt_tick:pkt_ylim);
yline(pkt12_lim, Color="#5e5e5e", LineStyle="--", Label=sprintf("%d Hz", pkt12_lim), LabelHorizontalAlignment="left");
legend(orbits_labels_sf12, Location="northeast");

subplot(3, 2, 4);
grid on; grid minor;
title(sprintf("Doppler rate, SX12xx, SF11, h = %d km", altitude));
ylabel("Doppler rate (Hz/T_{pkt})");
xlabel("Sim time");
xlim([dl_g dr_g]);
ylim([0, pkt_ylim]);
yticks(0:pkt_tick:pkt_ylim);
yline(pkt11_lim, Color="#5e5e5e", LineStyle="--", Label=sprintf("%d Hz", pkt11_lim), LabelHorizontalAlignment="left");
legend(orbits_labels_sf11, Location="northeast");

subplot(3, 2, 6);
grid on; grid minor;
title(sprintf("Doppler rate, SX12xx, SF10, h = %d km", altitude));
ylabel("Doppler rate (Hz/T_{pkt})");
xlabel("Sim time");
xlim([dl_g dr_g]);
ylim([0, pkt_ylim]);
yticks(0:pkt_tick:pkt_ylim);
yline(pkt10_lim, Color="#5e5e5e", LineStyle="--", Label=sprintf("%d Hz", pkt10_lim), LabelHorizontalAlignment="left");
legend(orbits_labels_sf10, Location="northeast");
