%% SATELLITE LINK BUDGET ANALYSIS
% Juan Del Pino Mena
% June 2025
%
% Link budget estimator for LEO satellites. Provides a complete link budget.
% This script computes access intervals in a given simulation time, azimuth, elevation,
% range, latency, Doppler frequency shift, FSPL losses, atmospheric losses, received power
% and CNR, across several orbits altitudes.
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
symbolperiod_s = 0.03277;  % [s] Symbol period
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
tx_loss_gs_dB = 2;  % [dB] TX system power loss in the ground station
rx_loss_gs_dB = 2;  % [dB] RX system power loss in the ground station


% ----------------------------------------------------------------------------------------
% Satellite

% Power, gain and losses

tx_power_sat_dBm = 33;  % [dBm] Transmission power in the satellite tx
tx_ant_gain_sat_dB = 2;  % [dBi] Antenna Gain in the sat
tx_loss_sat_dB = 2;  % [dB] System power loss in TX in the satellite
rx_loss_sat_dB = 2;  % [dB] System power loss in RX in the satellite


% ----------------------------------------------------------------------------------------
% Simulation time, stop time and step

dl = datetime('25-Sep-2024 02:25:00', 'Format', 'yyyy-MM-dd HH:mm:ss', TimeZone="UTC");
dr = dl + minutes(30); % datetime('25-Sep-2024 00:20:00', 'Format', 'yyyy-MM-dd HH:mm:ss', TimeZone="UTC");
sampleTime = 1;  % [s] Simulation step


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create scenario

scenario = satelliteScenario(dl, dr, sampleTime);  % Satellite scenario

altitudes = 300:50:450;
loop_iter = 0;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEGIN FOR LOOP
for altitude=altitudes

loop_iter = loop_iter + 1;


% ----------------------------------------------------------------------------------------
% Satellite

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
%% Compute access to the satellite from the ground station

ac_sat_gs = access(sat, gs);  % Compute access to the satellite
ac_sat_intervals = accessIntervals(ac_sat_gs);  % Access intervals from GS

[ac_sat_plot_data, ac_sat_time] = accessStatus(ac_sat_gs);  % Plot the access intervals
ac_sat_plot_data = double(ac_sat_plot_data);  % Casting for the operation below

ac_sat_plot_data_loop(loop_iter, :) = ac_sat_plot_data;


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


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Azimuth, elevation and range between satellite and ground station

% Compute azimuth, elevation, range
[azimuth_sat_gs_deg, elevation_sat_gs_deg, range_sat_gs_m] = aer(gs, sat);  % [deg,deg,m]

range_sat_gs_full_m = range_sat_gs_m;  % Save a copy of the complete vector for CNR calc.
aer_sat_gs_time = ac_sat_time;  % time vector for plotting

% Discard data when not in line of sight
azimuth_sat_gs_deg(~logical(ac_sat_plot_data)) = NaN;
elevation_sat_gs_deg(~logical(ac_sat_plot_data)) = NaN;
range_sat_gs_m(~logical(ac_sat_plot_data)) = NaN;

azimuth_sat_gs_deg_loop(loop_iter, :) = azimuth_sat_gs_deg;
elevation_sat_gs_deg_loop(loop_iter, :) = elevation_sat_gs_deg;
range_sat_gs_m_loop(loop_iter, :) = range_sat_gs_m;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Latency and Doppler Shift
% These functions are only available on MatLab R2023a onwards
% Dynamic - computes latency and doppler shift over the complete simulation time

[latency_delay, latency_time] = latency(sat, gs);
[doppler_fshift, doppler_time, doppler_info] = dopplershift(sat, gs, Frequency=freq_Hz);

latency_delay_loop(loop_iter, :) = latency_delay;
doppler_fshift_loop(loop_iter, :) = doppler_fshift;
doppler_time_loop(loop_iter, :) = doppler_time;
doppler_info_loop(loop_iter, :) = doppler_info;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% P618 propagation losses - atmospheric losses
% Estimate the atmospheric losses expected near Valencia.
% Static - does not recalculate for every simulation point

% ------------------------------------------------------------------------------
% Download and extract the digital maps, if not available on path

maps = exist("maps.mat","file");
p836 = exist("p836.mat","file");
p837 = exist("p837.mat","file");
p840 = exist("p840.mat","file");
matFiles = [maps p836 p837 p840];

if ~all(matFiles)
    if ~exist("ITURDigitalMaps.tar.gz","file")
        url = "https://www.mathworks.com/supportfiles/spc/P618/ITURDigitalMaps.tar.gz";
        websave("ITURDigitalMaps.tar.gz",url);
        untar("ITURDigitalMaps.tar.gz");
    else
        untar("ITURDigitalMaps.tar.gz");
    end
    addpath(cd)
end

% ------------------------------------------------------------------------------
% Configuration object

% KNOWN PROBLEMS:
% - The lowest frequency the model admits is 1 GHz.
% - Where is the polarization loss?
% - The cross-polarization discrimination prediction method is valid for the frequency
%   values in the range 4 to 55 GHz. For less than 4 GHz, the 4 GHz point will be used.

p618cfg = p618Config(...
    Frequency = 1e9, ...  % [Hz] Carrier frequency
    ElevationAngle = elevation_min_angle, ...  % [degrees] Elevation Angle, worst-case
    Latitude = gs_lat_N, ...  % Degrees North
    Longitude = gs_lon_E, ...  % Degrees East
    GasAnnualExceedance = 0.1, ...  % Average annual time percentage of excess for gas att
    CloudAnnualExceedance = 0.1, ...  % Avg annual time percentage of excess for cloud att
    RainAnnualExceedance = 0.001, ...  % Avg annual time percentage of excess for rain att
    ScintillationAnnualExceedance = 0.01, ...  % Avg annual time % of excess scintillation
    TotalAnnualExceedance = 0.001, ...  % Avg annual time % of excess for total att
    PolarizationTiltAngle = 0, ...  % [degrees, -90 to 90]  Polarization tilt angle
    AntennaDiameter = 1, ...  % [m] Physical diameter of the antenna. Default = 1 m
    AntennaEfficiency = 0.5 ...  % [0-1] Antenna efficiency
);

% ------------------------------------------------------------------------------
% Calculate Earth-space propagation losses, cross-polarization discrimination and sky
% noise temperature with the above configuration

[p618_atm_loss_dB, p618_xpol_discr_dB, p618_temp_sky_K] = p618PropagationLosses(p618cfg);

% p618_atm_loss_dB.Ag - Gaseous attenuation (dB)
% p618_atm_loss_dB.Ac - Cloud and fog attenuation (dB)
% p618_atm_loss_dB.Ar - Rain attenuation (dB)
% p618_atm_loss_dB.As - Attenuation due to tropospheric scintillation (dB)
% p618_atm_loss_dB.At - Total atmospheric attenuation (dB)
% p618_xpol_discr_dB  - Cross-polarization discrimination (dB) not exceeded for the
%                    percentage of the RainAnnualExceedance.
% p618_temp_sky_K     - Sky noise temperature (K) at the ground station antenna.

% ------------------------------------------------------------------------------
% Polarization loss

p618_pol_loss = 3;  % [dB] Polarization loss. Worst-case: 3 dB

% ------------------------------------------------------------------------------
% Total losses

p618_loss_total = p618_atm_loss_dB.At + p618_pol_loss;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Link budget analysis: Received Power

% EIRP: Equivalent Isotropic Radiated Power: EIRP = P_tx + G_ant_tx
% Pre_rx_loss: Total loss before the receiver input in the receiver system. Sum of feeder
%              loss, radome loss, loss due to polarization mismatch, etc.

% PISO: Received isotropic power in the antenna: PISO_rx = EIRP - FSPL
% PRI:  Received power at receiver input: PRI_rx = EIRP_tx - FSPL + G_ant_rx - Pre_rx_loss

% The PRI result is the most interesting. But since the PRI does not account for the p618
% losses, we will substract them from the received power

% ------------------------------------------------------------------------------
% Uplink

link_UL = link(gs_transmitter, sat_receiver);  % Create link
[link_UL_PISO_dBW, link_UL_PRI_dBW, link_UL_time] = sigstrength(link_UL);  % Signal power

% link_UL_rxpower_dBm = link_UL_PRI_dBW + 30 - p618_loss_total;  % [dBm] UL received power
link_UL_rxpower_dBm = link_UL_PRI_dBW + 30 - p618_loss_total + tx_ant_gain_sat_dB + tx_ant_gain_gs_dB;  % [dBm] UL received power

% ------------------------------------------------------------------------------
% Downlink

link_DL = link(sat_transmitter, gs_receiver);  % Create link
[link_DL_PISO_dBW, link_DL_PRI_dBW, link_DL_time] = sigstrength(link_DL);  % Signal power

% link_DL_rxpower_dBm = link_DL_PRI_dBW + 30 - p618_loss_total;  % [dBm] DL received power
link_DL_rxpower_dBm = link_DL_PRI_dBW + 30 - p618_loss_total + tx_ant_gain_sat_dB + tx_ant_gain_gs_dB;  % [dBm] DL received power

link_UL_rxpower_dBm_loop(loop_iter, :) = link_UL_rxpower_dBm;
link_DL_rxpower_dBm_loop(loop_iter, :) = link_DL_rxpower_dBm;


% PFD

pfd_loop(loop_iter, :) = (tx_power_sat_dBm - 30) + (10 * log10(4000 / bandwidth_Hz)) + ...
                         (10 * log10(1 ./ (4 * pi * range_sat_gs_m.^2))) - p618_loss_total;


end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END FOR LOOP
%% GRAPHS

close all;

dl_g = dl + minutes(9);  % limits
dr_g = dr - minutes(9);

orbits_labels = {};
i = 1;
for altitude=altitudes
    orbits_labels{i} = sprintf("%d km", altitude);
    i = i + 1;
end


% ------------------------------------------------------------------------------
% Plot points in time where there is access to the satellite.

figure();
plot(ac_sat_time, ac_sat_plot_data_loop); hold on;
title("Observation intervals: line-of-sight with $\varepsilon >= \varepsilon_{min}$" + ...
      ". [1: access; 0: no access]", interpreter="latex");
xlabel("Sim time");
ylabel("Access status");

ylim([-0.5, 1.5]);
xlim([dl_g dr_g]);
grid on; grid minor;

legend(orbits_labels, Location="south");


% ------------------------------------------------------------------------------
% Plot azimuth, elevation, range

figure();

% azimuth_sat_gs_deg_loop(azimuth_sat_gs_deg_loop <= 180) = 180 + azimuth_sat_gs_deg_loop(azimuth_sat_gs_deg_loop <= 180);

% subplot(3, 1, 1);
% plot(aer_sat_gs_time, azimuth_sat_gs_deg_loop); hold on;
% grid on; grid minor;
% xlabel("Simulation time");
% ylabel("Azimuth (degrees)");
% title("Ground Station Azimuth");
% xlim([dl_g dr_g]);
% ylim([220, 430]);
% yticks(220:30:430);
% legend(orbits_labels, Location="northeast");

subplot(1, 2, 1);
plot(aer_sat_gs_time, elevation_sat_gs_deg_loop); hold on;
grid on; grid minor;
xlabel("Sim time");
ylabel("Elevation (degrees)");
title("Ground Station Elevation");
xlim([dl_g dr_g]);
ylim([0, 90]);
yticks(0:15:90);
legend(orbits_labels, Location="northeast");

subplot(1, 2, 2);
plot(aer_sat_gs_time, range_sat_gs_m_loop * 1e-3); hold on;
grid on; grid minor;
xlabel("Sim time");
ylabel("Distance (km)");
title("Range between GS and Sat");
xlim([dl_g dr_g]);
ylim([0, 2000]);
yticks(0:500:2000);


% ------------------------------------------------------------------------------
% Plot latency, doppler shift

figure();

subplot(1, 2, 1);
plot(latency_time, latency_delay_loop .* 1e3); hold on;  % in milliseconds
xlim([latency_time(1), latency_time(end)]);
title("Radio link one-way latency");
xlabel("Sim time");
ylabel("Latency (ms)");
grid on; grid minor;
ylim([0, 7]);
xlim([dl_g dr_g]);
yticks(0:1:7);
% legend(orbits_labels, Location="southwest");

subplot(1, 2, 2);
plot(doppler_time, doppler_fshift_loop * 1e-3); hold on;  % kHz
xlim([doppler_time(1), doppler_time(end)]);
title("Radio link Doppler shift");
xlabel("Sim time");
ylabel("Doppler Shift (kHz)");
grid on; grid minor;
ylim([-25, 25]);
xlim([dl_g dr_g]);
yticks(-25:5:25);
legend(orbits_labels, Location="northeast");


% ------------------------------------------------------------------------------
%% Plot received power

figure();

ylim_p = [-142,-112];
ylim_p_tick = 3;

subplot(1, 2, 1);
plot(link_DL_time, link_DL_rxpower_dBm_loop);
hold on;
grid on; grid minor;
xlabel("Sim time");
ylabel("Received power (dBm)");
title("Received signal strength, downlink");
yline(rx_limit_sensitivity_gs, Color="#5e5e5e", LineStyle="--", Label="S SF12"); hold on;
yline(rx_limit_sensitivity_gs+3, Color="#5e5e5e", LineStyle="--", Label="S SF11"); hold on;
yline(rx_limit_sensitivity_gs+6, Color="#5e5e5e", LineStyle="--", Label="S SF10"); hold on;
ylim(ylim_p);
xlim([dl_g dr_g]);
yticks(min(ylim_p):ylim_p_tick:max(ylim_p));

legend(orbits_labels, Location="south");

subplot(1, 2, 2);

plot(link_UL_time, zeros(1, length(link_UL_time)), '-', Color="#000000"); hold on;
plot(link_UL_time, zeros(1, length(link_UL_time)), '-.',  Color="#000000");

plot(link_UL_time, link_UL_rxpower_dBm_loop(1, :) + 8, LineStyle="-", Color="#0072bd"); % 22 dBm
plot(link_UL_time, link_UL_rxpower_dBm_loop(2, :) + 8, LineStyle="-", Color="#d95319");
plot(link_UL_time, link_UL_rxpower_dBm_loop(3, :) + 8, LineStyle="-", Color="#edb120");
plot(link_UL_time, link_UL_rxpower_dBm_loop(4, :) + 8, LineStyle="-", Color="#7e2f8e");

plot(link_UL_time, link_UL_rxpower_dBm_loop(1, :), LineStyle="-.", Color="#0072bd"); % 14 dBm
plot(link_UL_time, link_UL_rxpower_dBm_loop(2, :), LineStyle="-.", Color="#d95319");
plot(link_UL_time, link_UL_rxpower_dBm_loop(3, :), LineStyle="-.", Color="#edb120");
plot(link_UL_time, link_UL_rxpower_dBm_loop(4, :), LineStyle="-.", Color="#7e2f8e");

grid on; grid minor;
xlabel("Sim time");
ylabel("Received power (dBm)");
title("Received signal strength, uplink");
yline(rx_limit_sensitivity_sat, Color="#5e5e5e", LineStyle="--", Label="S SF12"); hold on;
yline(rx_limit_sensitivity_sat+3, Color="#5e5e5e", LineStyle="--", Label="S SF11"); hold on;
yline(rx_limit_sensitivity_sat+6, Color="#5e5e5e", LineStyle="--", Label="S SF10"); hold on;
ylim(ylim_p);
xlim([dl_g dr_g]);
yticks(min(ylim_p):ylim_p_tick:max(ylim_p));

power_labels = {'22 dBm', '14 dBm'};
legend(power_labels, Location="northeast");


% ------------------------------------------------------------------------------
%% Plot power to respect PFD

figure();

ylim_p = [-156,-136];
ylim_p_tick = 3;

plot(ac_sat_time, pfd_loop);
hold on;
grid on; grid minor;
xlabel("Sim time");
ylabel("Power Flux Density (dB(W/m^2·4 kHz))");
title("Power Flux Density on Earth's surface (P_{sat} = 33 dBm)");
ylim(ylim_p);
xlim([dl_g dr_g]);
yticks(min(ylim_p):ylim_p_tick:max(ylim_p));
yline(-142, Color="#5e5e5e", LineStyle="--", Label="-142 dB(W/m^2·4 kHz)"); hold on;

legend(orbits_labels, Location="south");
