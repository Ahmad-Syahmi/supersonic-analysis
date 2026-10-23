clear; clc; close all;

%% Geometry, calibration, and channel mapping
x_taps = [24.7, 36.6, 48.3, 61.0, 73.7, 86.4, 99.1, 118.0, 209.5, 300.1];
x_pitot = [24.7, 36.6, 48.3, 61.0, 73.7, 86.4, 99.1];
A_over_As = [1.06, 1.00, 1.05, 1.15, 1.23, 1.27, 1.28, 1.30, 1.30, 1.30];
x_throat = 36.6;
gamma = 1.4;
channel_map = [1, 2, 3, 5, 4, 6, 7, 8, 9, 11, 10, 12];
[~, sort_order] = sort(channel_map);
K_gain = 1.7838e4;
u_K_gain = 0.5;                % Pa/V, from +/-0.00005e4 Pa/V
p_atm = 100050;               % Pa, 100.05 kPa measured atmospheric pressure
u_p_atm = 410;                % Pa, 0.41 kPa standard uncertainty
transducer_accuracy = 0.0025;  % fraction of measured gauge pressure

data_dir = fileparts(mfilename('fullpath'));
sup_files = dir(fullfile(data_dir, 'supersonic_*_data_*.dat'));
sub_files = dir(fullfile(data_dir, 'subsonic_*_data_*.dat'));
sup_offsets = dir(fullfile(data_dir, 'supersonic_*_offset_*.dat'));
sub_offsets = dir(fullfile(data_dir, 'subsonic_*_offset_*.dat'));

if isempty(sup_files) || isempty(sub_files)
    error('Both subsonic and supersonic *_data_*.dat files are required.');
end

sup = process_files(sup_files, sup_offsets, data_dir, sort_order, ...
    K_gain, u_K_gain, p_atm, u_p_atm, transducer_accuracy, gamma);
sub = process_files(sub_files, sub_offsets, data_dir, sort_order, ...
    K_gain, u_K_gain, p_atm, u_p_atm, transducer_accuracy, gamma);
if numel(sup) ~= numel(x_pitot) || numel(sub) ~= numel(x_pitot)
    error('Expected exactly seven Port01-Port07 data files per regime.');
end

fprintf('Supersonic: %d files x %d static taps = %d plotted points.\n', ...
    numel(sup), numel(x_taps), numel(sup) * numel(x_taps));
fprintf('Subsonic: %d files x %d static taps = %d plotted points.\n', ...
    numel(sub), numel(x_taps), numel(sub) * numel(x_taps));
fprintf('Static taps plotted: 1-10. Channel 12 is the total-pressure reference.\n');
fprintf('Mach points plotted at Pitot leading edges: %d per flow condition.\n', ...
    numel(x_pitot));

%% Theoretical isentropic solution
x_dense = linspace(min(x_taps), max(x_taps), 400);
AR_dense = interp1(x_taps, A_over_As, x_dense, 'pchip');
M_theory_sub = zeros(size(x_dense));
M_theory_sup = zeros(size(x_dense));
idx_th = find(x_dense >= x_throat, 1);
for i = 1:numel(x_dense)
    AR = AR_dense(i);
    area_fun = @(M) (1 / M) * ((2 / (gamma + 1)) * ...
        (1 + (gamma - 1) / 2 * M^2))^((gamma + 1) / (2 * (gamma - 1))) - AR;
    M_theory_sub(i) = fzero(area_fun, [1e-4, 1.0]);
    if i < idx_th
        M_theory_sup(i) = fzero(area_fun, [1e-4, 1.0]);
    else
        M_theory_sup(i) = fzero(area_fun, [1.0, 3.0]);
    end
end
p_rat_theory_sub = (1 + (gamma - 1) / 2 .* M_theory_sub.^2).^(-gamma / (gamma - 1));
p_rat_theory_sup = (1 + (gamma - 1) / 2 .* M_theory_sup.^2).^(-gamma / (gamma - 1));
p0_reservoir = mean([sub.p0]);
p_theory_sub = p_rat_theory_sub * p0_reservoir;
p_theory_sup = p_rat_theory_sup * p0_reservoir;
p0_pitot_theory_sub = repmat(p0_reservoir, size(x_dense));
rayleigh_ratio = rayleigh_pitot_ratio(M_theory_sup, gamma);
p0_pitot_theory_sup = p_theory_sup .* rayleigh_ratio;
p0_pitot_theory_sup(x_dense < x_throat) = p0_reservoir;

%% Plot the five requested distributions
set(groot, 'defaultAxesFontName', 'Times New Roman', ...
    'defaultTextFontName', 'Times New Roman');
static_taps = 1:numel(x_taps);
sup_static = vertcat(sup.p_static);
sub_static = vertcat(sub.p_static);
sup_static_u = vertcat(sup.u_p_static);
sub_static_u = vertcat(sub.u_p_static);
sup_static_mean = mean(sup_static, 1);
sub_static_mean = mean(sub_static, 1);
sup_static_u_mean = sqrt(sum(sup_static_u.^2, 1)) / numel(sup);
sub_static_u_mean = sqrt(sum(sub_static_u.^2, 1)) / numel(sub);

figure('Name', 'Static pressure over atmospheric pressure', ...
    'Color', [1 1 1]);
plot(x_dense, p_rat_theory_sup, '-', 'Color', [0.85 0.1 0.1], ...
    'LineWidth', 1.5, 'DisplayName', 'Supersonic isentropic theory');
hold on;
plot(x_dense, p_rat_theory_sub, '--', 'Color', [0.1 0.3 0.8], ...
    'LineWidth', 1.3, 'DisplayName', 'Subsonic isentropic theory');
errorbar(x_taps(static_taps), sup_static_mean / p_atm, ...
    sup_static_u_mean / p_atm, '-s', ...
    'Color', [0.85 0.1 0.1], 'MarkerFaceColor', [0.85 0.1 0.1], ...
    'DisplayName', 'Supersonic (10 static taps)');
errorbar(x_taps(static_taps), sub_static_mean / p_atm, ...
    sub_static_u_mean / p_atm, '-o', ...
    'Color', [0.1 0.3 0.8], 'MarkerFaceColor', [0.1 0.3 0.8], ...
    'DisplayName', 'Subsonic (10 static taps)');
xline(x_throat, 'k:', 'LineWidth', 1.2, 'DisplayName', 'Nozzle throat');
xlabel('Axial Position, x [mm]');
ylabel('Normalized static pressure, p / p_{atm}');
title('Static Pressure Distribution, p / p_{atm}');
grid on; box on; xlim([20 310]); ylim([0.2 1.05]);
legend('Location', 'best');

plot_pressure_ratio(x_pitot, sub, p_rat_theory_sub, x_dense, ...
    x_throat, 'Subsonic', [0.1 0.3 0.8], 'o');
plot_pressure_ratio(x_pitot, sup, p_rat_theory_sup, x_dense, ...
    x_throat, 'Supersonic', [0.85 0.1 0.1], 's');
plot_mach(x_pitot, sub, M_theory_sub, x_dense, ...
    x_throat, 'Subsonic', [0.1 0.3 0.8], 'o');
plot_mach(x_pitot, sup, M_theory_sup, x_dense, ...
    x_throat, 'Supersonic', [0.85 0.1 0.1], 's');
plot_stagnation_pressure(x_pitot, sup, sub, p0_pitot_theory_sup, ...
    p0_pitot_theory_sub, x_dense, x_throat, p_atm);

function plot_pressure_ratio(x, results, theory_ratio, theory_x, ...
    x_throat, regime, color, marker)
port_numbers = 1:numel(results);
measured_ratio = arrayfun(@(result, port) result.p_ratio(port), ...
    results, port_numbers);
measured_u = arrayfun(@(result, port) result.u_ratio(port), ...
    results, port_numbers);

figure('Name', sprintf('%s p/p0', regime), 'Color', [1 1 1]);
plot(theory_x, theory_ratio, '-', 'Color', color, 'LineWidth', 1.5, ...
    'DisplayName', 'Isentropic theory');
hold on;
errorbar(x, measured_ratio, measured_u, ['-' marker], ...
    'Color', color, 'MarkerFaceColor', color, 'LineWidth', 1.1, ...
    'DisplayName', sprintf('%s measured (7 points)', regime));
xline(x_throat, 'k:', 'LineWidth', 1.2, 'DisplayName', 'Nozzle throat');
xlabel('Axial Position, x [mm]');
ylabel('Static pressure ratio, p / p_0');
title(sprintf('%s Static Pressure Ratio, p / p_0', regime));
grid on; box on; xlim([min(x) - 5, max(x) + 5]); ylim([0 1.05]);
legend('Location', 'best');
end

function plot_mach(x, results, theory_mach, theory_x, ...
    x_throat, regime, color, marker)
measured_mach = [results.mach_pitot];
measured_u = [results.u_mach_pitot];

figure('Name', sprintf('%s Mach number', regime), 'Color', [1 1 1]);
plot(theory_x, theory_mach, '-', 'Color', color, 'LineWidth', 1.5, ...
    'DisplayName', 'Isentropic theory');
hold on;
errorbar(x, measured_mach, measured_u, ['-' marker], ...
    'Color', color, 'MarkerFaceColor', color, 'LineWidth', 1.1, ...
    'DisplayName', sprintf('%s measured (7 points)', regime));
xline(x_throat, 'k:', 'LineWidth', 1.2, 'DisplayName', 'Nozzle throat');
xlabel('Axial Position, x [mm]');
ylabel('Mach number, M');
title(sprintf('%s Mach Number Distribution', regime));
grid on; box on; xlim([min(x) - 5, max(x) + 5]); ylim([0 1.6]);
legend('Location', 'best');
end

function plot_stagnation_pressure(x, sup, sub, theory_sup, theory_sub, ...
    theory_x, x_throat, p_atm)
figure('Name', 'Stagnation pressure', 'Color', [1 1 1]);
plot(theory_x, theory_sup / p_atm, '--', 'Color', [0.85 0.1 0.1], ...
    'LineWidth', 1.4, ...
    'DisplayName', 'Supersonic post-shock Pitot theory');
hold on;
plot(theory_x, theory_sub / p_atm, ':', 'Color', [0.1 0.3 0.8], ...
    'LineWidth', 1.4, ...
    'DisplayName', 'Subsonic reservoir theory');
errorbar(x, [sup.p0] / p_atm, [sup.u_p0] / p_atm, '-s', ...
    'Color', [0.85 0.1 0.1], 'MarkerFaceColor', [0.85 0.1 0.1], ...
    'LineWidth', 1.1, ...
    'DisplayName', 'Measured supersonic Pitot stagnation pressure');
errorbar(x, [sub.p0] / p_atm, [sub.u_p0] / p_atm, '-o', ...
    'Color', [0.1 0.3 0.8], 'MarkerFaceColor', [0.1 0.3 0.8], ...
    'LineWidth', 1.1, ...
    'DisplayName', 'Measured subsonic stagnation pressure');
xline(x_throat, 'k:', 'LineWidth', 1.2, 'DisplayName', 'Nozzle throat');
xlabel('Pitot axial position, x [mm]');
ylabel('Normalized stagnation pressure, p_0 / p_{atm}');
title('Stagnation Pressure Distribution');
grid on; box on;
xlim([min(x) - 5, max(x) + 5]);
legend('Location', 'best');
end

function results = process_files(data_files, offset_files, data_dir, ...
    sort_order, K_gain, u_K_gain, p_atm, u_p_atm, ...
    transducer_accuracy, gamma)
offset_data = [];
offset_sample_intervals = zeros(1, numel(offset_files));
for k = 1:numel(offset_files)
    offset_raw = load_dat_file(fullfile(data_dir, offset_files(k).name));
    offset_sample_intervals(k) = median(diff(offset_raw(:, 1)));
    offset_data = [offset_data; offset_raw]; %#ok<AGROW>
end
if isempty(offset_data)
    error('No offset file was found for this regime.');
end
offset_mean = mean(offset_data(:, 2:end), 1);
if any(~isfinite(offset_sample_intervals)) || ...
        any(offset_sample_intervals <= 0)
    error('Invalid time spacing in an offset data file.');
end
offset_sample_interval = median(offset_sample_intervals);
offset_u_voltage_95 = zeros(1, 12);
for channel = 1:12
    [~, offset_n_independent] = calculate_integral_time_scale( ...
        offset_data(:, channel + 1), offset_sample_interval);
    offset_u_voltage_95(channel) = 1.96 * ...
        std(offset_data(:, channel + 1), 0) / sqrt(offset_n_independent);
end
offset_u_pressure_95 = K_gain * offset_u_voltage_95;
offset_u_pressure_95 = offset_u_pressure_95(sort_order);
results = struct('name', {}, 'p_ratio', {}, 'u_ratio', {}, ...
    'mach_pitot', {}, 'u_mach_pitot', {}, 'mach_method', {}, ...
    'port_number', {}, ...
    'p0', {}, 'u_p0', {}, ...
    'p_static', {}, 'u_p_static', {}, ...
    'integral_time_scale', {}, 'n_independent', {}, ...
    'u_p_combined_95', {});

for k = 1:numel(data_files)
    raw = load_dat_file(fullfile(data_dir, data_files(k).name));
    samples = raw(:, 2:end);
    n = size(samples, 1);
    sample_interval = median(diff(raw(:, 1)));
    if ~isfinite(sample_interval) || sample_interval <= 0
        error('Invalid time spacing in %s.', data_files(k).name);
    end
    p_gauge = (mean(samples, 1) - offset_mean) * K_gain;
    p_sd = std(samples, 0, 1) * K_gain;
    p_gain_u = abs(mean(samples, 1) - offset_mean) * u_K_gain;
    p_gauge = p_gauge(sort_order);
    p_sd = p_sd(sort_order);
    p_gain_u = p_gain_u(sort_order);

    p_abs = p_atm + p_gauge;
    p_accuracy_u = transducer_accuracy * abs(p_gauge);
    p_precision_95 = zeros(1, 12);
    integral_time_scale = zeros(1, 12);
    n_independent = zeros(1, 12);
    for channel = 1:12
        [integral_time_scale(channel), n_independent(channel)] = ...
            calculate_integral_time_scale( ...
            samples(:, sort_order(channel)), sample_interval);
        p_precision_95(channel) = 1.96 * p_sd(channel) / ...
            sqrt(n_independent(channel));
    end
    u_accuracy_precision_95 = hypot(p_accuracy_u, p_precision_95);
    u_p = hypot(u_accuracy_precision_95, ...
        hypot(p_gain_u, hypot(offset_u_pressure_95, u_p_atm)));
    p0 = p_abs(12);
    u_p0 = u_p(12);
    if p0 <= 0
        error('Non-positive total pressure in %s.', data_files(k).name);
    end

    ratio_all = p_abs ./ p0;
    u_ratio_all = ratio_all .* hypot(u_p ./ p_abs, u_p0 / p0);
    results(k).name = data_files(k).name;
    results(k).p0 = p0;
    results(k).u_p0 = u_p0;
    results(k).integral_time_scale = integral_time_scale;
    results(k).n_independent = n_independent;
    results(k).u_p_combined_95 = u_p;
    results(k).p_static = p_abs(1:10);
    results(k).u_p_static = u_p(1:10);
    port_number = parse_port_number(data_files(k).name);
    if port_number < 1 || port_number > 7
        error('Port number %d in %s is outside the expected range 1-7.', ...
            port_number, data_files(k).name);
    end
    local_tap_index = port_number;
    if contains(data_files(k).name, 'supersonic_')
        [results(k).mach_pitot, results(k).u_mach_pitot, ...
            results(k).mach_method] = ...
            calc_M_rayleigh_pitot(p_abs(local_tap_index), p0, ...
            u_p(local_tap_index), u_p0, gamma);
        mach = results(k).mach_pitot;
        ratio_all(local_tap_index) = ...
            (1 + (gamma - 1) / 2 * mach^2)^(-gamma / (gamma - 1));
        d_ratio_d_mach = -gamma * mach * ...
            (1 + (gamma - 1) / 2 * mach^2)^ ...
            (-(2 * gamma - 1) / (gamma - 1));
        u_ratio_all(local_tap_index) = ...
            abs(d_ratio_d_mach) * results(k).u_mach_pitot;
    else
        results(k).mach_pitot = calc_M_isen(p_abs(local_tap_index), p0, gamma);
        results(k).u_mach_pitot = calc_M_isen_uncertainty( ...
            p_abs(local_tap_index), p0, u_p(local_tap_index), u_p0, gamma);
        results(k).mach_method = 'isentropic';
    end
    results(k).p_ratio = ratio_all(1:10);
    results(k).u_ratio = u_ratio_all(1:10);
    results(k).port_number = port_number;
    fprintf('  %s: %d samples, 10 static points\n', data_files(k).name, n);
end

[~, port_order] = sort([results.port_number]);
results = results(port_order);
end

function data = load_dat_file(file_path)
fid = fopen(file_path, 'r');
if fid == -1
    error('Could not open file: %s', file_path);
end
cleanup = onCleanup(@() fclose(fid));
line = fgetl(fid);
while ischar(line)
    if startsWith(strtrim(line), 'Time')
        break;
    end
    line = fgetl(fid);
end
if ~ischar(line)
    error('Could not find the Time header in %s.', file_path);
end
data = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f', ...
    'Delimiter', '\t', 'CollectOutput', true);
data = data{1};
if isempty(data) || size(data, 2) ~= 13
    error('Expected a Time column and 12 DAQ channels in %s.', file_path);
end
end

function mach = calc_M_isen(static_pressure, stagnation_pressure, gamma)
pressure_term = (stagnation_pressure ./ static_pressure).^((gamma - 1) / gamma) - 1;
mach = sqrt(max(0, (2 / (gamma - 1)) .* pressure_term));
end

function u_mach = calc_M_isen_uncertainty(static_pressure, stagnation_pressure, ...
    u_static, u_stagnation, gamma)
mach = calc_M_isen(static_pressure, stagnation_pressure, gamma);
if mach <= 0 || static_pressure <= 0 || stagnation_pressure <= 0
    u_mach = 0;
    return;
end
exponent = (gamma - 1) / gamma;
d_mach_dp = -stagnation_pressure^exponent * ...
    static_pressure^(-exponent - 1) / (gamma * mach);
d_mach_dp0 = stagnation_pressure^(exponent - 1) * ...
    static_pressure^(-exponent) / (gamma * mach);
u_mach = hypot(d_mach_dp * u_static, d_mach_dp0 * u_stagnation);
end

function [mach, u_mach, method] = calc_M_rayleigh_pitot(static_pressure, pitot_pressure, ...
    u_static, u_pitot, gamma)
if static_pressure <= 0 || pitot_pressure <= 0
    error('Rayleigh-Pitot pressures must be positive.');
end
ratio = pitot_pressure / static_pressure;
if ratio < 1
    error('Rayleigh-Pitot pressure ratio must be at least one; got %.6g.', ratio);
end
sonic_ratio = rayleigh_pitot_ratio(1, gamma);
if ratio < sonic_ratio
    mach = calc_M_isen(static_pressure, pitot_pressure, gamma);
    u_mach = calc_M_isen_uncertainty(static_pressure, pitot_pressure, ...
        u_static, u_pitot, gamma);
    method = 'isentropic fallback (R < R(M=1))';
    return;
end
mach = solve_rayleigh_pitot(ratio, gamma);
u_ratio = ratio * hypot(u_pitot / pitot_pressure, u_static / static_pressure);
d_ratio_d_mach = rayleigh_pitot_ratio(mach, gamma) * ...
    (2 * gamma / ((gamma - 1) * mach) - ...
    4 * gamma * mach / ((gamma - 1) * ...
    (2 * gamma * mach^2 - (gamma - 1))));
u_mach = u_ratio / abs(d_ratio_d_mach);
method = 'Rayleigh-Pitot';
end

function mach = solve_rayleigh_pitot(ratio, gamma)
upper = 3;
while rayleigh_pitot_ratio(upper, gamma) < ratio && upper < 20
    upper = upper * 2;
end
if rayleigh_pitot_ratio(upper, gamma) < ratio
    error('Rayleigh-Pitot ratio %.6g exceeds solver range.', ratio);
end
mach = fzero(@(candidate) rayleigh_pitot_ratio(candidate, gamma) - ratio, ...
    [1, upper]);
end

function ratio = rayleigh_pitot_ratio(mach, gamma)
denominator = 2 * gamma * mach.^2 - (gamma - 1);
ratio = (((gamma + 1) / 2) .* mach.^2).^(gamma / (gamma - 1)) .* ...
    ((gamma + 1) ./ denominator).^(1 / (gamma - 1));
end

function port_number = parse_port_number(file_name)
token = regexp(file_name, 'Port(\d+)_data', 'tokens', 'once');
if isempty(token)
    error('Could not extract a port number from filename: %s', file_name);
end
port_number = str2double(token{1});
end

function [integral_time_scale, n_independent] = ...
    calculate_integral_time_scale(signal, sample_interval)
signal = signal - mean(signal);
n_samples = numel(signal);
variance = sum(signal.^2);
if variance == 0 || n_samples < 2
    integral_time_scale = 0;
    n_independent = n_samples;
    return;
end

autocorrelation = zeros(n_samples, 1);
for lag = 0:n_samples - 1
    autocorrelation(lag + 1) = ...
        sum(signal(1:n_samples - lag) .* signal(1 + lag:n_samples)) / ...
        variance;
end

first_nonpositive = ...
    find(autocorrelation(2:end) <= 0, 1, 'first') + 1;
if isempty(first_nonpositive)
    last_lag = n_samples;
else
    last_lag = first_nonpositive;
end

lags = (0:last_lag - 1)' * sample_interval;
integral_time_scale = trapz(lags, autocorrelation(1:last_lag));
integral_time_scale = max(0, integral_time_scale);
n_independent = max(1, floor(n_samples * sample_interval / ...
    max(2 * integral_time_scale, sample_interval)));
end
