clear; clc;

data_dir = fileparts(mfilename('fullpath'));
gamma = 1.4;
K_gain = 1.7838e4;
u_K_gain = 0.5;                % Pa/V, from +/-0.00005e4 Pa/V
p_atm = 101325;
u_p_atm = 150;
transducer_accuracy = 0.0025;  % fraction of measured pressure magnitude

x_taps = [24.7, 36.6, 48.3, 61.0, 73.7, 86.4, 99.1, 118.0, ...
    209.5, 300.1];
channel_map = [1, 2, 3, 5, 4, 6, 7, 8, 9, 11, 10, 12];
[~, sort_order] = sort(channel_map);

regimes = ["subsonic", "supersonic"];
output = table();

for regime = regimes
    data_files = dir(fullfile(data_dir, regime + "_*_data_*.dat"));
    offset_files = dir(fullfile(data_dir, regime + "_*_offset_*.dat"));
    if isempty(data_files)
        warning('No %s data files found.', regime);
        continue;
    end

    offset_data = [];
    for k = 1:numel(offset_files)
        offset_data = [offset_data; load_dat_file(fullfile( ...
            data_dir, offset_files(k).name))];
    end
    if isempty(offset_data)
        error('No valid offset data found for %s.', regime);
    end
    offset_mean = mean(offset_data(:, 2:end), 1);
    offset_sample_interval = median(diff(offset_data(:, 1)));
    if ~isfinite(offset_sample_interval) || offset_sample_interval <= 0
        error('Invalid time spacing in the %s offset data.', regime);
    end
    offset_u_voltage = zeros(1, 12);
    for channel = 1:12
        [~, offset_n_independent] = calculate_integral_time_scale( ...
            offset_data(:, channel + 1), offset_sample_interval);
        offset_u_voltage(channel) = std(offset_data(:, channel + 1), 0) / ...
            sqrt(offset_n_independent);
    end
    offset_u_pressure = K_gain * offset_u_voltage;
    offset_u_pressure = offset_u_pressure(sort_order);

    for k = 1:numel(data_files)
        file_name = data_files(k).name;
        port_number = parse_port_number(file_name);
        if port_number < 1 || port_number > 7
            error('Port number %d in %s is outside the expected range 1-7.', ...
                port_number, file_name);
        end
        raw = load_dat_file(fullfile(data_dir, file_name));
        samples = raw(:, 2:end);
        n_samples = size(samples, 1);

        voltage_mean = mean(samples, 1);
        voltage_sd = std(samples, 0, 1);
        voltage_sem = voltage_sd / sqrt(n_samples);
        time = raw(:, 1);
        sample_interval = median(diff(time));
        if ~isfinite(sample_interval) || sample_interval <= 0
            error('Invalid time spacing in %s.', file_name);
        end

        p_gauge = (voltage_mean - offset_mean) * K_gain;
        p_sd = voltage_sd * K_gain;
        p_sem = voltage_sem * K_gain;
        p_gain_u = abs(voltage_mean - offset_mean) * u_K_gain;

        p_gauge = p_gauge(sort_order);
        p_sd = p_sd(sort_order);
        p_sem = p_sem(sort_order);
        p_gain_u = p_gain_u(sort_order);

        p_abs = p_atm + p_gauge;
        if any(p_abs <= 0)
            error('Non-positive absolute pressure in %s.', file_name);
        end
        u_p_sd = hypot(p_sd, hypot(offset_u_pressure, u_p_atm));
        u_p_sem = hypot(p_sem, hypot(offset_u_pressure, u_p_atm));
        p_accuracy_u = transducer_accuracy * abs(p_gauge);
        p_precision_95 = zeros(1, 12);
        n_independent = zeros(1, 12);
        integral_time_scale = zeros(1, 12);
        for channel = 1:12
            [integral_time_scale(channel), n_independent(channel)] = ...
                calculate_integral_time_scale( ...
                samples(:, sort_order(channel)), sample_interval);
            p_precision_95(channel) = 1.96 * p_sd(channel) / ...
                sqrt(n_independent(channel));
        end
        u_accuracy_precision_95 = hypot(p_accuracy_u, p_precision_95);
        u_p_combined_95 = hypot( ...
            u_accuracy_precision_95, ...
            hypot(p_gain_u, hypot(1.96 * offset_u_pressure, u_p_atm)));

        p0 = p_abs(12);
        u_p0_sd = u_p_sd(12);
        u_p0_sem = u_p_sem(12);
        if p0 <= 0
            error('Non-positive port-12 absolute pressure in %s.', file_name);
        end

        p_ratio = p_abs ./ p0;
        u_ratio_sd = p_ratio .* hypot(u_p_sd ./ p_abs, u_p0_sd / p0);
        u_ratio_sem = p_ratio .* hypot(u_p_sem ./ p_abs, u_p0_sem / p0);

        if regime == "supersonic"

            mach = NaN(12, 1);
            u_mach_sd = NaN(12, 1);
            u_mach_sem = NaN(12, 1);
            local_tap_index = port_number;
            [mach(local_tap_index), mach_method] = calc_mach_rayleigh_pitot( ...
                p_abs(local_tap_index), p0, gamma);
            u_mach_sd(local_tap_index) = calc_mach_rayleigh_uncertainty( ...
                p_abs(local_tap_index), p0, ...
                u_p_sd(local_tap_index), u_p0_sd, gamma);
            u_mach_sem(local_tap_index) = calc_mach_rayleigh_uncertainty( ...
                p_abs(local_tap_index), p0, ...
                u_p_sem(local_tap_index), u_p0_sem, gamma);
        else
            mach = calc_mach(p_abs, p0, gamma);
            u_mach_sd = calc_mach_uncertainty(p_abs, p0, ...
                u_p_sd, u_p0_sd, gamma);
            u_mach_sem = calc_mach_uncertainty(p_abs, p0, ...
                u_p_sem, u_p0_sem, gamma);
            mach_method = "isentropic";
        end

        [~, stem] = fileparts(file_name);
        tap_number = (1:12)';
        x_position = [x_taps, NaN, NaN]';
        rows = table( ...
            repmat(regime, 12, 1), repmat(string(stem), 12, 1), ...
            repmat(n_samples, 12, 1), repmat(port_number, 12, 1), ...
            tap_number, x_position, ...
            p_gauge', p_sd', p_sem', p_abs', u_p_sd', u_p_sem', ...
            offset_u_pressure', ...
            integral_time_scale', n_independent', p_precision_95', ...
            p_accuracy_u', p_gain_u', u_accuracy_precision_95', ...
            u_p_combined_95', ...
            p_ratio', u_ratio_sd', u_ratio_sem', mach', ...
            u_mach_sd', u_mach_sem', repmat(mach_method, 12, 1), ...
            'VariableNames', {'Regime', 'DataFile', 'Samples', 'Port', 'Tap', ...
            'X_mm', 'GaugePressure_Pa', 'PressureSD_Pa', 'PressureSEM_Pa', ...
            'AbsolutePressure_Pa', 'AbsolutePressureU_SD_Pa', ...
            'AbsolutePressureU_SEM_Pa', ...
            'OffsetMeanU_Pa', ...
            'IntegralTimeScale_s', 'N_independent', 'Precision95_Pa', ...
            'TransducerAccuracy_Pa', 'GainUncertainty_Pa', ...
            'AccuracyPlusPrecisionU95_Pa', 'CombinedU95_Pa', ...
            'PressureRatio', 'PressureRatioU_SD', ...
            'PressureRatioU_SEM', ...
            'Mach', 'MachU_SD', 'MachU_SEM', 'MachMethod'});
        output = [output; rows];
    end
end

output_file = fullfile(data_dir, 'pressure_tap_uncertainty_by_file.csv');
writetable(output, output_file);
fprintf('Wrote %d rows to %s\n', height(output), output_file);

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

function mach = calc_mach(static_pressure, stagnation_pressure, gamma)
pressure_term = (stagnation_pressure ./ static_pressure).^((gamma - 1) / gamma) - 1;
mach = sqrt(max(0, (2 / (gamma - 1)) .* pressure_term));
end

function u_mach = calc_mach_uncertainty(static_pressure, stagnation_pressure, ...
    u_static, u_stagnation, gamma)
mach = calc_mach(static_pressure, stagnation_pressure, gamma);
exponent = (gamma - 1) / gamma;
valid = mach > 0 & static_pressure > 0 & stagnation_pressure > 0;
d_mach_dp = zeros(size(mach));
d_mach_dp0 = zeros(size(mach));
d_mach_dp(valid) = -stagnation_pressure.^exponent .* ...
    static_pressure(valid).^(-exponent - 1) ./ (gamma * mach(valid));
d_mach_dp0(valid) = stagnation_pressure.^(exponent - 1) .* ...
    static_pressure(valid).^(-exponent) ./ (gamma * mach(valid));
u_mach = hypot(d_mach_dp .* u_static, d_mach_dp0 .* u_stagnation);
end

function [mach, method] = calc_mach_rayleigh_pitot(static_pressure, pitot_pressure, gamma)
if any(static_pressure <= 0) || pitot_pressure <= 0
    error('Rayleigh-Pitot pressures must be positive.');
end
ratio = pitot_pressure ./ static_pressure;
if any(ratio < 1)
    error('Rayleigh-Pitot pressure ratios must be at least one.');
end
sonic_ratio = rayleigh_pitot_ratio(1, gamma);
mach = zeros(size(ratio));
subsonic = ratio < sonic_ratio;
mach(subsonic) = calc_mach(static_pressure(subsonic), pitot_pressure, gamma);
mach(~subsonic) = arrayfun(@(value) solve_rayleigh_pitot(value, gamma), ...
    ratio(~subsonic));
method = 'Rayleigh-Pitot';
if any(subsonic)
    method = 'isentropic fallback (R < R(M=1))';
end
end

function u_mach = calc_mach_rayleigh_uncertainty(static_pressure, pitot_pressure, ...
    u_static, u_pitot, gamma)
mach = calc_mach_rayleigh_pitot(static_pressure, pitot_pressure, gamma);
ratio = pitot_pressure ./ static_pressure;
u_ratio = ratio .* hypot(u_pitot ./ pitot_pressure, ...
    u_static ./ static_pressure);
sonic_ratio = rayleigh_pitot_ratio(1, gamma);
subsonic = ratio < sonic_ratio;
u_mach = zeros(size(ratio));
u_mach(subsonic) = calc_mach_uncertainty( ...
    static_pressure(subsonic), pitot_pressure, ...
    u_static(subsonic), u_pitot, gamma);
valid = ~subsonic;
d_ratio_d_mach = rayleigh_pitot_ratio(mach, gamma) .* ...
    (2 * gamma ./ ((gamma - 1) .* mach) - ...
    4 * gamma .* mach ./ ((gamma - 1) .* ...
    (2 * gamma .* mach.^2 - (gamma - 1))));
u_mach(valid) = u_ratio(valid) ./ abs(d_ratio_d_mach(valid));
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

% Normalized autocorrelation for non-negative lags. The integral ends at
% the first non-positive correlation, which avoids integrating noise-driven
% negative lobes at long lags.
autocorrelation = zeros(n_samples, 1);
for lag = 0:n_samples - 1
    autocorrelation(lag + 1) = ...
        sum(signal(1:n_samples - lag) .* signal(1 + lag:n_samples)) / ...
        variance;
end
first_nonpositive = find(autocorrelation(2:end) <= 0, 1, 'first') + 1;
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
