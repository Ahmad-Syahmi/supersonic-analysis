clear; clc; close all;

%% Geometry, calibration, and channel mapping
x_taps = [24.7, 36.6, 48.3, 61.0, 73.7, 86.4, 99.1, 118.0, 209.5, 300.1];
x_pitot = [36.6, 48.3, 61.0, 73.7, 86.4, 99.1, 118.0];
A_over_As = [1.06, 1.00, 1.05, 1.15, 1.23, 1.27, 1.28, 1.30, 1.30, 1.30];
x_throat = 36.6;
gamma = 1.4;
channel_map = [1, 2, 3, 5, 4, 6, 7, 8, 9, 11, 10, 12];
[~, sort_order] = sort(channel_map);
K_gain = 1.7838e4;
p_atm = 101325;
u_p_atm = 150;

data_dir = fileparts(mfilename('fullpath'));
sup_files = dir(fullfile(data_dir, 'supersonic_*_data_*.dat'));
sub_files = dir(fullfile(data_dir, 'subsonic_*_data_*.dat'));
sup_offsets = dir(fullfile(data_dir, 'supersonic_*_offset_*.dat'));
sub_offsets = dir(fullfile(data_dir, 'subsonic_*_offset_*.dat'));

if isempty(sup_files) || isempty(sub_files)
    error('Both subsonic and supersonic *_data_*.dat files are required.');
end

sup = process_files(sup_files, sup_offsets, data_dir, sort_order, ...
    K_gain, p_atm, u_p_atm, gamma);
sub = process_files(sub_files, sub_offsets, data_dir, sort_order, ...
    K_gain, p_atm, u_p_atm, gamma);
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

%% Plot pressure ratio and Mach number
set(groot, 'defaultAxesFontName', 'Times New Roman', ...
    'defaultTextFontName', 'Times New Roman');
figure('Name', 'Pressure and Mach distributions', ...
    'Units', 'inches', 'Position', [1, 1, 10.5, 7.0], ...
    'Color', [1 1 1]);
tab_group = uitabgroup();
pressure_tab = uitab(tab_group, 'Title', 'Static pressure');
mach_tab = uitab(tab_group, 'Title', 'Mach number');
stagnation_tab = uitab(tab_group, 'Title', 'Stagnation pressure');

axes('Parent', pressure_tab);
theory_sup_patm = p_rat_theory_sup * mean([sup.p0]) / p_atm;
theory_sub_patm = p_rat_theory_sub * mean([sub.p0]) / p_atm;
plot(x_dense, theory_sup_patm, 'r-', 'LineWidth', 1.5, ...
    'DisplayName', 'Theory: Supersonic'); hold on;
plot(x_dense, theory_sub_patm, 'b--', 'LineWidth', 1.3, ...
    'DisplayName', 'Theory: Subsonic');
sup_colors = lines(max(numel(sup), numel(sub)));
selected_ports = [3, 5, 7];
for k = selected_ports
    errorbar(x_taps, sup(k).p_static / p_atm, ...
        sup(k).u_p_static / p_atm, '-s', ...
        'Color', sup_colors(k, :), 'MarkerFaceColor', sup_colors(k, :), ...
        'LineWidth', 1.0, 'CapSize', 3, ...
        'DisplayName', sprintf('Supersonic Port%02d', k));
end
for k = selected_ports
    errorbar(x_taps, sub(k).p_static / p_atm, ...
        sub(k).u_p_static / p_atm, '-o', ...
        'Color', sup_colors(k, :), 'MarkerFaceColor', sup_colors(k, :), ...
        'LineWidth', 1.0, 'CapSize', 3, ...
        'DisplayName', sprintf('Subsonic Port%02d', k));
end
xline(x_throat, 'k:', 'LineWidth', 1.2, 'DisplayName', 'Nozzle throat');
xlabel('Axial Position, x [mm]');
ylabel('Normalized static pressure, p / p_{atm}');
title('(a) Static Pressure Distribution, p / p_{atm}');
xlim([20 310]); ylim([0.2 1.05]); legend('Location', 'eastoutside', 'FontSize', 7);

axes('Parent', mach_tab);
% Supersonic points use the Rayleigh-Pitot shock correction; subsonic
% points use the isentropic total-to-static relation.
plot(x_dense, M_theory_sup, 'r-', 'LineWidth', 1.5, ...
    'DisplayName', 'Theory: Supersonic'); hold on;
plot(x_dense, M_theory_sub, 'b--', 'LineWidth', 1.3, ...
    'DisplayName', 'Theory: Subsonic');
mach_ports = [1:5, 7];
errorbar(x_pitot(mach_ports), [sup(mach_ports).mach_pitot], ...
    [sup(mach_ports).u_mach_pitot], '-s', ...
    'Color', [0.85 0.1 0.1], ...
    'MarkerFaceColor', [0.85 0.1 0.1], 'LineWidth', 1.1, ...
    'DisplayName', 'Supersonic (Rayleigh-Pitot, Port06 omitted)');
errorbar(x_pitot(mach_ports), [sub(mach_ports).mach_pitot], ...
    [sub(mach_ports).u_mach_pitot], '-o', ...
    'Color', [0.1 0.3 0.8], ...
    'MarkerFaceColor', [0.1 0.3 0.8], 'LineWidth', 1.1, ...
    'DisplayName', 'Subsonic (isentropic, Port06 omitted)');
xline(x_throat, 'k:', 'LineWidth', 1.2, 'DisplayName', 'Nozzle throat');
yline(1.0, 'k--', 'Alpha', 0.5, 'HandleVisibility', 'off');
xlabel('Axial Position, x [mm]'); ylabel('Mach number, M');
title('(b) Mach Number Distribution'); grid on; box on;
xlim([20 310]); ylim([0 1.6]); legend('Location', 'eastoutside', 'FontSize', 7);

axes('Parent', stagnation_tab);
errorbar(x_pitot, [sup.p0] / p_atm, [sup.u_p0] / p_atm, '-s', ...
    'Color', [0.85 0.1 0.1], ...
    'MarkerFaceColor', [0.85 0.1 0.1], 'LineWidth', 1.1, ...
    'DisplayName', 'Supersonic');
hold on;
errorbar(x_pitot, [sub.p0] / p_atm, [sub.u_p0] / p_atm, '-o', ...
    'Color', [0.1 0.3 0.8], ...
    'MarkerFaceColor', [0.1 0.3 0.8], 'LineWidth', 1.1, ...
    'DisplayName', 'Subsonic');
xlabel('Pitot axial position, x [mm]');
ylabel('Normalized stagnation pressure, p_0 / p_{atm}');
title('Stagnation Pressure Distribution');
grid on; box on;
xlim([min(x_pitot) - 5, max(x_pitot) + 5]);
legend('Location', 'best');

function results = process_files(data_files, offset_files, data_dir, ...
    sort_order, K_gain, p_atm, u_p_atm, gamma)
offset_data = [];
for k = 1:numel(offset_files)
    offset_data = [offset_data; load_dat_file(fullfile( ...
        data_dir, offset_files(k).name))]; %#ok<AGROW>
end
if isempty(offset_data)
    error('No offset file was found for this regime.');
end
offset_mean = mean(offset_data(:, 2:end), 1);
results = struct('name', {}, 'p_ratio', {}, 'u_ratio', {}, ...
    'mach_pitot', {}, 'u_mach_pitot', {}, 'mach_method', {}, ...
    'port_number', {}, ...
    'p0', {}, 'u_p0', {}, ...
    'p_static', {}, 'u_p_static', {});

for k = 1:numel(data_files)
    raw = load_dat_file(fullfile(data_dir, data_files(k).name));
    samples = raw(:, 2:end);
    n = size(samples, 1);
    p_gauge = (mean(samples, 1) - offset_mean) * K_gain;
    p_sd = std(samples, 0, 1) * K_gain;
    p_gauge = p_gauge(sort_order);
    p_sd = p_sd(sort_order);

    p_abs = p_atm + p_gauge;
    u_p = hypot(p_sd, u_p_atm);
    p0 = p_abs(12);
    u_p0 = u_p(12);
    if p0 <= 0
        error('Non-positive total pressure in %s.', data_files(k).name);
    end

    ratio_all = p_abs ./ p0;
    u_ratio_all = ratio_all .* hypot(u_p ./ p_abs, u_p0 / p0);
    results(k).name = data_files(k).name;
    results(k).p_ratio = ratio_all(1:10);
    results(k).u_ratio = u_ratio_all(1:10);
    results(k).p0 = p0;
    results(k).u_p0 = u_p0;
    results(k).p_static = p_abs(1:10);
    results(k).u_p_static = u_p(1:10);
    port_number = parse_port_number(data_files(k).name);
    if port_number < 1 || port_number > 7
        error('Port number %d in %s is outside the expected range 1-7.', ...
            port_number, data_files(k).name);
    end
    local_tap_index = port_number + 1;
    if contains(data_files(k).name, 'supersonic_')
        [results(k).mach_pitot, results(k).u_mach_pitot, ...
            results(k).mach_method] = ...
            calc_M_rayleigh_pitot(p_abs(local_tap_index), p0, ...
            u_p(local_tap_index), u_p0, gamma);
    else
        results(k).mach_pitot = calc_M_isen(p_abs(local_tap_index), p0, gamma);
        results(k).u_mach_pitot = calc_M_isen_uncertainty( ...
            p_abs(local_tap_index), p0, u_p(local_tap_index), u_p0, gamma);
        results(k).mach_method = 'isentropic';
    end
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
