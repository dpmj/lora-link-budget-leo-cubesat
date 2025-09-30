%% SATELLITE LINK BUDGET ANALYSIS
% Juan Del Pino Mena
% June 2025
%
% Link budget estimator for LEO satellites, sweeping the orbit altitude.
% Variant for establishing a Power Flux Density limit.
% This script only calculates the received power and the carrier-to-noise ratio.
%
% This script requires MatLab >= R2023a and the satellite communications toolbox.


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
symbolrate_Sps = 1/symbolperiod_s;  % [Sps] Symbol rate
bitrate_bps = 293;  % [b/s] Transmission bitrate

% ----------------------------------------------------------------------------------------
% Transceiver limits (only used for plots)

rx_limit_sensitivity_gs = -137;  % [dBm] The minimum admitted signal power in the gs
rx_limit_cnr_gs = -20;  % [dB] The minimum SNR / CNR the gs transceiver admits

rx_limit_sensitivity_sat = -141;  % [dBm] The minimum admitted signal power in the sat
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
% Warning: This script will plot results for tx_power_gs_dBm = 14 dBm and
% tx_power_gs_dBm + 8 = 22 dBm, the maximum power allowed for the SRD (EU868) band.
% See "Plots" section.

tx_ant_gain_gs_dB = 2;  % [dBi] Antenna Gain in the GS
tx_loss_gs_dB = 2;  % [dB] TX system power loss in the ground station
rx_loss_gs_dB = 2;  % [dB] RX system power loss in the ground station

% Noise

rx_nf_gs_dB = 6;  % [dB] System noise figure in RX in the ground station

ant_ambient_temp_gs_K = 300;  % [K] Antenna ambient temperature in the ground station
ant_noise_temp_gs_K = 300;  % [K] Antenna noise temperature in the ground station

% [dB/K] G/T in the gs rx
rx_GNTR_gs_dB_K = tx_ant_gain_gs_dB - rx_nf_gs_dB - 10 * log10(ant_ambient_temp_gs_K + ...
    (ant_noise_temp_gs_K - ant_ambient_temp_gs_K) * 10^(-0.1 * rx_nf_gs_dB));  % [dB/K]


% ----------------------------------------------------------------------------------------
% Satellite

% Power, gain and losses

tx_power_sat_dBm = 33;  % [dBm] Transmission power in the satellite tx
% Warning: This script will reduce tx_power_sat_dBm according to the orbit altitude as
% specified in the variable tx_power_sat_pfd_reduction to comply with the Power Flux
% Density limits on Earth.

tx_ant_gain_sat_dB = 2;  % [dBi] Antenna Gain in the sat
tx_loss_sat_dB = 2;  % [dB] System power loss in TX in the satellite
rx_loss_sat_dB = 2;  % [dB] System power loss in RX in the satellite

% Noise

rx_nf_sat_dB = 3;  % [dB] System noise figure in RX in the satellite

ant_ambient_temp_sat_K = 350;  % [K] Antenna ambient temperature in the satellite
ant_noise_temp_sat_K = 350;  % [K] Antenna noise temperature in the satellite

% [dB/K] G/T in the sat rx
rx_GNTR_sat_dB_K = tx_ant_gain_sat_dB - rx_nf_sat_dB ...
    - 10 * log10(ant_ambient_temp_sat_K ...
    + (ant_noise_temp_sat_K - ant_ambient_temp_sat_K) * 10^(-0.1 * rx_nf_sat_dB));


% ----------------------------------------------------------------------------------------
% Simulation time, stop time and step

dl = datetime('25-Sep-2024 02:25:00', 'Format', 'yyyy-MM-dd HH:mm:ss', TimeZone="UTC");
dr = dl + minutes(30);
sampleTime = 1;  % [s] Simulation step


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scenario

scenario = satelliteScenario(dl, dr, sampleTime);  % Satellite scenario

altitudes = 300:50:450;
loop_iter = 0;
%                                 300,     350,     400,     450      [km]
tx_power_sat_pfd_reduction = [-4.4129, -3.0572, -1.8847, -0.8518];  % [dB]
% How much to reduce the satellite transmission power to comply with PFD regulations


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEGIN MAIN FOR LOOP
for altitude=altitudes

loop_iter = loop_iter + 1;


% ----------------------------------------------------------------------------------------
% Satellite Orbit

semiMajorAxis = (R_E + altitude) * 1e3;  % [m]
eccentricity = 0;
inclination = 50;  % [deg]
rightAscensionOfAscendingNode = 0;  % [deg]
argumentOfPeriapsis = 0;  % [deg]
trueAnomaly = 0;  % [deg]

sat = satellite(scenario, ...
                semiMajorAxis, ...  % [m]
                eccentricity, ...
                inclination, ...  % [deg]
                rightAscensionOfAscendingNode,...  % [deg]
                argumentOfPeriapsis, ...  % [deg]
                trueAnomaly, ...  % [deg]
                Name="sat");  % Satellite name


% ------------------------------------------------------------------------------
% Ground Station

gs = groundStation(scenario, Name=gs_name, Latitude=gs_lat_N, Longitude=gs_lon_E, ...
                   Altitude=gs_altitude, MinElevationAngle=elevation_min_angle);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute access to the satellite from the ground station

ac_sat_gs = access(sat, gs);  % Compute access to the satellite
ac_sat_intervals = accessIntervals(ac_sat_gs);  % Access intervals from GS
[ac_sat_plot_data, ac_sat_time] = accessStatus(ac_sat_gs);  % Plot the access intervals


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

% Caviats:
% - This configuration will vary according to the region of the ground station
% - The lowest frequency the model admits is 1 GHz.
% - The polarization loss is manually picked as 3 dB
% - The cross-polarization discrimination prediction method is valid for the frequency
%   values in the range 4 to 55 GHz. For less than 4 GHz, the 4 GHz point is used.

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
%                       percentage of the RainAnnualExceedance.
% p618_temp_sky_K     - Sky noise temperature (K) at the ground station antenna.

% ------------------------------------------------------------------------------
% Polarization loss

p618_pol_loss = 3;  % [dB] Polarization loss. Worst-case: 3 dB (loss of one component)

% ------------------------------------------------------------------------------
% Total losses

p618_loss_total = p618_atm_loss_dB.At + p618_pol_loss;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CNR: Carrier-to-noise ratio for configured satellite link budget parameters
% Computes CNR in every simulation point

N_points_sim = length(range_sat_gs_m);

ul_CNR = zeros(N_points_sim, 1);
dl_CNR = zeros(N_points_sim, 1);

for i = 1:1:N_points_sim

    % UPLINK: Ground Station -> Satellite
    ul_CNR_cfg = satelliteCNRConfig(...
        TransmitterPower = tx_power_gs_dBm - 30, ...  % [dBW] Transmit power
        TransmitterSystemLoss = tx_loss_gs_dB, ...  % [dB] Losses in the transmitter
        TransmitterAntennaGain = tx_ant_gain_gs_dB, ...  % [dBi] Antenna gain
        Distance = range_sat_gs_m(i) * 1e-3, ...  % [km] Distance. Varies.
        Frequency = freq_Hz * 1e-9, ...  % [GHz] Center frequency
        MiscellaneousLoss = p618_loss_total, ...  % [dB] Misc attenuation
        GainToNoiseTemperatureRatio = rx_GNTR_sat_dB_K, ...  % [dB/K] Receiver G/T ratio
        ReceiverSystemLoss = rx_loss_sat_dB, ...  % [dB] Receiver system loss
        BitRate = bitrate_bps * 1e-6, ...  % [Mbps] Bit rate
        SymbolRate = symbolrate_Sps * 1e-6, ...  % [MSps] Symbol rate
        Bandwidth = bandwidth_Hz * 1e-6 ...  % [MHz] Bandwidth
    );
    [ul_CNR(i), ~] = satelliteCNR(ul_CNR_cfg);

    % DOWNLINK: Satellite -> Ground Station
    dl_CNR_cfg = satelliteCNRConfig(...
        TransmitterPower = tx_power_sat_dBm - 30, ...  % [dBW] Transmit power
        TransmitterSystemLoss = tx_loss_sat_dB, ...  % [dB] Losses in the transmitter
        TransmitterAntennaGain = tx_ant_gain_sat_dB, ...  % [dBi] Antenna gain
        Distance = range_sat_gs_m(i) * 1e-3, ...  % [km] Distance. Varies.
        Frequency = freq_Hz * 1e-9, ...  % [GHz] Center frequency
        MiscellaneousLoss = -tx_power_sat_pfd_reduction(loop_iter) + p618_loss_total, ...  % [dB] Misc attenuation
        GainToNoiseTemperatureRatio = rx_GNTR_gs_dB_K, ...  % [dB/K] Receiver G/T ratio
        ReceiverSystemLoss = rx_loss_gs_dB, ...  % [dB] Receiver system loss
        BitRate = bitrate_bps * 1e-6, ...  % [Mbps] Bit rate
        SymbolRate = symbolrate_Sps * 1e-6, ...  % [MSps] Symbol rate
        Bandwidth = bandwidth_Hz * 1e-6 ...  % [MHz] Bandwidth
    );
    [dl_CNR(i), ~] = satelliteCNR(dl_CNR_cfg);

end

% Filter data out of the line of sight
ul_CNR(~logical(ac_sat_plot_data)) = NaN;
dl_CNR(~logical(ac_sat_plot_data)) = NaN;

cnr_time = ac_sat_time;  % time vector for plotting

ul_CNR_loop(loop_iter, :) = ul_CNR;  % Cannot allocate due to type not being double
dl_CNR_loop(loop_iter, :) = dl_CNR;


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
link_DL_rxpower_dBm = link_DL_PRI_dBW + 30 - p618_loss_total + tx_ant_gain_sat_dB + tx_ant_gain_gs_dB + tx_power_sat_pfd_reduction(loop_iter);  % [dBm] DL received power

link_UL_rxpower_dBm_loop(loop_iter, :) = link_UL_rxpower_dBm;
link_DL_rxpower_dBm_loop(loop_iter, :) = link_DL_rxpower_dBm;


% PFD

pfd_loop(loop_iter, :) = (tx_power_sat_dBm - 30) + (10 * log10(4000 / bandwidth_Hz)) + ...
                         (10 * log10(1 ./ (4 * pi * range_sat_gs_m.^2))) - p618_loss_total + tx_power_sat_pfd_reduction(loop_iter);


end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END MAIN FOR LOOP
%% GRAPHS

% Prepare the plot data

dl_g = dl + minutes(9);  % visualization limits
dr_g = dr - minutes(9);

orbits_labels = {};  % label for plots
i = 1;
for altitude=altitudes
    orbits_labels{i} = sprintf("%d km, %.1f dBm", altitude, tx_power_sat_dBm + tx_power_sat_pfd_reduction(i));
    i = i + 1;
end


%% Plots

% ------------------------------------------------------------------------------
% Plot Carrier to Noise Ratio

figure();

ylim_p = [-30,0];  % plot y axis limit
ylim_p_tick = 3;  % plot y axis tick step

subplot(1, 2, 1);
plot(cnr_time, dl_CNR_loop);
hold on;
grid on; grid minor;
xlabel("Sim time");
ylabel("SNR (dB)");
title("Signal-to-Noise Ratio, downlink");
yline(rx_limit_cnr_gs, Color="#5e5e5e", LineStyle="--", Label="SNR SF12"); hold on;
yline(rx_limit_cnr_gs+2.5, Color="#5e5e5e", LineStyle="--", Label="SNR SF11"); hold on;
yline(rx_limit_cnr_gs+5, Color="#5e5e5e", LineStyle="--", Label="SNR SF10"); hold on;
ylim(ylim_p);
xlim([dl_g dr_g]);
yticks(min(ylim_p):ylim_p_tick:max(ylim_p));
leg = legend(orbits_labels, Location="southwest");
title(leg, 'Altitude, Sat TX power')

subplot(1, 2, 2);
plot(cnr_time, -100*ones(1, length(cnr_time)), '-', Color="#000000"); hold on;
plot(cnr_time, -100*ones(1, length(cnr_time)), '-.',  Color="#000000");

plot(cnr_time, ul_CNR_loop(1, :), LineStyle="-.", Color="#0072bd"); % 14 dBm
plot(cnr_time, ul_CNR_loop(2, :), LineStyle="-.", Color="#d95319");
plot(cnr_time, ul_CNR_loop(3, :), LineStyle="-.", Color="#edb120");
plot(cnr_time, ul_CNR_loop(4, :), LineStyle="-.", Color="#7e2f8e");

plot(cnr_time, ul_CNR_loop(1, :) + 8, LineStyle="-", Color="#0072bd"); % 22 dBm
plot(cnr_time, ul_CNR_loop(2, :) + 8, LineStyle="-", Color="#d95319");
plot(cnr_time, ul_CNR_loop(3, :) + 8, LineStyle="-", Color="#edb120");
plot(cnr_time, ul_CNR_loop(4, :) + 8, LineStyle="-", Color="#7e2f8e");

grid on; grid minor;
xlabel("Sim time");
ylabel("SNR (dB)");
title("Signal-to-Noise Ratio, uplink");
yline(rx_limit_cnr_sat, Color="#5e5e5e", LineStyle="--", Label="SNR SF12"); hold on;
yline(rx_limit_cnr_sat+2.5, Color="#5e5e5e", LineStyle="--", Label="SNR SF11"); hold on;
yline(rx_limit_cnr_sat+5, Color="#5e5e5e", LineStyle="--", Label="SNR SF10"); hold on;
ylim(ylim_p);
xlim([dl_g dr_g]);
yticks(min(ylim_p):ylim_p_tick:max(ylim_p));

power_labels = {'22 dBm', '14 dBm'};
leg = legend(power_labels, Location="northeast");
title(leg, 'GS TX power')


% ------------------------------------------------------------------------------
% Plot received power

figure();

ylim_p = [-146,-116];
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

leg = legend(orbits_labels, Location="southwest");
title(leg, 'Altitude, Sat TX power')

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
leg = legend(power_labels, Location="northeast");
title(leg, 'GS TX power')


% ------------------------------------------------------------------------------
% Plot PFD

figure();

ylim_p = [-156,-136];
ylim_p_tick = 3;

plot(ac_sat_time, pfd_loop);
hold on;
grid on; grid minor;
xlabel("Sim time");
ylabel("PFD (dB(W/m^2·4 kHz))");
title("Power Flux Density on Earth's surface");
ylim(ylim_p);
xlim([dl_g dr_g]);
yticks(min(ylim_p):ylim_p_tick:max(ylim_p));
yline(-142, Color="#5e5e5e", LineStyle="--", Label="Limit = -142 dB(W/m^2·4 kHz)", ...
    LabelHorizontalAlignment='left'); hold on;

leg = legend(orbits_labels, Location="northeast");
title(leg, 'Altitude, Sat TX power')



