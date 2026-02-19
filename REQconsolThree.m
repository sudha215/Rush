clearvars
close all
clc

% paths to responce data (to be replaced by wav imports)
addpath('C:\Users\FinnAlderson\Matlab Cambridge\Text exports');

% filename = 'lyngdorf p1.txt';
filename = 'LRX new stands.txt';
data = readmatrix(filename, 'FileType', 'text', 'NumHeaderLines', 14);
fs = 44100;
freq = round(data(:,1), 5, "decimals");
mag  = round(data(:,2), 5, "decimals");
max_freq = 2000;
last_idx = find(freq <= max_freq, 1, 'last');
freq = freq(1:last_idx);
mag  = mag(1:last_idx);

% ensure column vectors
freq = freq(:);
mag  = mag(:);

clear last_idx

%% Key Parameters

minus = 12;                 % -x dB point for Fc finding
mean_grad = -5.8;             % set gradient for target slope
LF_flat_freq = 60;          % target curve lower flat section start
% maxIters =  15;              % number of filter itterations

% Limiter perameters
f_ref = 20;                 % reference low frequency where limiter is 0 dB
f_cap = 250;                % frequency where you want the limiter to reach limit_max
limit_max = 6;              % maximum allowed boost (dB)


%% Fast smoothing

S_order = [1 6 36 24 3];                         % meaning 1/N octave smoothing where N is S_order
nOrders = length(S_order);
SM = zeros(length(freq), nOrders);
log2f = log2(freq);                              % precompute log2(freq)

for iO = 1:nOrders
    N = S_order(iO);
    fprintf('Processing 1/%d octave smoothing ...\n', N);

    width = 1/(2*N);                             % half-width in octaves for the Gaussian
    dlog = log2f - log2f.';                      % NxN matrix: element (i,j) = log2(freq_i) - log2(freq_j)
    weights = exp(-(dlog.^2) / (2*width^2));     % gaussian weights in log2-domain
    wsum = sum(weights,1);
    SM(:,iO) = (weights * mag) ./ wsum.';
end

avgSmooth = mean(SM, 2);
mag = SM(:,3);
log10f = log10(freq);
Npts = length(freq);

% preallocate
mag_varSmooth = zeros(Npts, 1);

% Compute variable smoothing
N = zeros(size(freq));
for i = 1:Npts
    x = log10f(i);

    if x <= 2
        N(i) = 48; % 1/48 octave below 100 Hz
    elseif x >= 4
        N(i) = 3; % 1/3 octave above 10 kHz
    elseif x <= 3
        % between 100 and 1k Hz: interpolate 1/48 → 1/6 octave
        val = 1/48 + (x - 2) * ((1/6 - 1/48) / (3 - 2));
        N(i) = 1 / val;
    else
        % between 1k and 10k Hz: interpolate 1/6 → 1/3 octave
        val = 1/6 + (x - 3) * ((1/3 - 1/6) / (4 - 3));
        N(i) = 1 / val;
    end
end

for i = 1:Npts
    width = 1 / (2 * N(i)); % Gaussian half-width in octaves
    dlog = log2f - log2f(i);
    w = exp(-(dlog.^2) / (2 * width^2));
    w = w / sum(w);
    mag_varSmooth(i) = sum(w .* mag);

end

figure('Name','Variable smoothing')
semilogx(freq,mag_varSmooth,'DisplayName','Variable smoothing','LineWidth',2)
hold on
semilogx(freq,SM(:,1),'DisplayName',sprintf('1/%d smoothing',S_order(5)))
semilogx(freq,SM(:,5),'DisplayName',sprintf('1/%d smoothing',S_order(5)))
semilogx(freq,mag,'DisplayName','import data')
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
grid on; legend show
xlim([20 2000])

%% find Fc

% =========================================================================  minus = 12;
smoothed = mag_varSmooth; % Use your variable smoothing curve

% Compute gradient in dB per decade
grad_log = gradient(smoothed) ./ gradient(log10(freq));
presmooth = grad_log;

% Smooth the gradient output
smooth_order = 12;
width = 1 / (2 * smooth_order);
grad_log_smooth = zeros(size(freq));

for i = 1:length(freq)
    weights = exp(-(log2(freq / freq(i))).^2 / (2 * width^2));
    grad_log_smooth(i) = sum(weights .* grad_log) / sum(weights);
end

grad_log = grad_log_smooth;

% Limit search range
freqRange = [30 200];
idxRange = (freq >= freqRange(1) & freq <= freqRange(2));

freq_sub = freq(idxRange);
grad_sub = grad_log(idxRange);

crossIdx = find(diff(sign(grad_sub))); % finds all places where the gradient crosses zero

% Compute zero-crossing frequencies
f_cross_all = zeros(size(crossIdx));
for i = 1:length(crossIdx)
    f1 = freq_sub(crossIdx(i));
    f2 = freq_sub(crossIdx(i) + 1);
    g1 = grad_sub(crossIdx(i));
    g2 = grad_sub(crossIdx(i) + 1);
    f_cross_all(i) = 10^(interp1([g1 g2], log10([f1 f2]), 0));
end


idx_mean_start = find(freq >= 400, 1, 'first');
idx_mean_end = find(freq >= 1000, 1, 'first');

mag_at_cross = interp1(freq, smoothed, f_cross_all, 'linear');

% added -10 to capture low magnitide, low order roll offs
threshold_dB = mean(avgSmooth(idx_mean_start:idx_mean_end)) - 10 ;

fprintf('Fc threshold is %.2f dB \n', threshold_dB);

valid_idx = mag_at_cross > threshold_dB;

% Keep only valid crossings
f_cross_valid = f_cross_all(valid_idx);
mag_cross_valid = mag_at_cross(valid_idx);

if isempty(f_cross_valid)
    warning('No valid zero-crossings found above %.1f dB threshold.', threshold_dB);
    f_zeroGrad = NaN;
else
    f_zeroGrad = f_cross_valid(1);
end

if ~isnan(f_zeroGrad)
    target_dB = interp1(freq, smoothed, f_zeroGrad, 'linear') - minus;
else
    target_dB = NaN;
end

% Find where the smoothed diff curve crosses that target line

idx_range = (freq >= 0 & freq <= 150);
freq_sub = freq(idx_range);
mag_sub  = smoothed(idx_range);
crossings = find(diff(sign(mag_sub - target_dB)));

f_cross = [];
for i = 1:length(crossings)
    f1 = freq_sub(crossings(i));
    f2 = freq_sub(crossings(i)+1);
    y1 = mag_sub(crossings(i));
    y2 = mag_sub(crossings(i)+1);
    f_cross(i) = f1 + (f2 - f1) * (target_dB - y1) / (y2 - y1);
end

figure('Name','Find Fhp','NumberTitle','off');
subplot(2,1,1);
semilogx(freq, grad_log, 'LineWidth', 1.5, 'DisplayName', 'Smoothed Gradient');
hold on;
semilogx(freq, presmooth, '--', 'Color', [0.6 0.6 0.6], 'DisplayName', 'Unsmoothed Gradient');
xlabel('Frequency (Hz)'); ylabel('Gradient (dB per decade)');
grid on;
yline(0, '--k', 'DisplayName', 'Zero line');
xlim([20 500]);
title('Smoothed Gradient and Zero Crossings');

if exist('f_cross_all','var')
    for i = 1:length(f_cross_all)
        color = [0.7 0.7 0.7];
        if valid_idx(i), color = [0.2 0.7 0.2]; end
        xline(f_cross_all(i), '--', 'Color', color, ...
            'DisplayName', sprintf('Cross %.1f Hz', f_cross_all(i)));
    end
end
if ~isnan(f_zeroGrad)
    xline(f_zeroGrad, 'r-', 'LineWidth', 1.2, 'DisplayName', sprintf('Selected Fc = %.1f Hz', f_zeroGrad));
end
legend show;

subplot(2,1,2);
semilogx(freq, smoothed, 'LineWidth', 1.5, 'DisplayName', 'Smoothed');
hold on;

% plot horizontal target line
hLine = yline(target_dB, '--', 'Color', [0.5 0.5 0.5], ...
    'DisplayName', sprintf('-%.fdB reference', minus));
yline(threshold_dB,'--','Color',[0.5 0.5 0.5],'DisplayName','Zero crossing threshold');


% Plot Fc point on the smoothed curve (y = smoothed at Fc)
if ~isnan(f_zeroGrad)
    mag_at_Fc = interp1(freq, smoothed, f_zeroGrad, 'pchip');
    semilogx(f_zeroGrad, mag_at_Fc, 'rx', 'MarkerSize',10,'MarkerFaceColor','r', ...
        'DisplayName',sprintf('Fc = %.1f Hz (on curve)', f_zeroGrad));
end

% Plot the -minus dB crossing points (frequencies where smoothed == target_dB)
if exist('f_cross','var') && ~isempty(f_cross)

    % highlight the first crossing
    semilogx(f_cross(1), target_dB, 'ro', 'MarkerSize',4, ...
        'DisplayName', sprintf('Selected -%.fdB at %.1f Hz', minus, f_cross(1)));
end

grid on; legend show;
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Smoothed Curve & Target Reference');
xlim([20 500]);


%% ------------------ Find dB grad for line section -----------------------

flat_start = 200;
% flat_start = f_cross(1,1)+150;
idx_flat_start = find(freq >= flat_start, 1, 'first');
gfs = crossings(1,1);

gf_freqs = freq(gfs:idx_flat_start);
gf_mag_int = SM(:,1);
gf_mag = gf_mag_int(gfs:idx_flat_start);


% gradient in dB per octave
grad_dB_per_oct = gradient(gf_mag) ./ gradient(log2(gf_freqs));
% mean_grad = mean(grad_dB_per_oct, 'omitnan');


% % set grad based on Dynaudios
% ========================================================================  mean_grad = -3;
fprintf('mean grad is: %.1f dB/oct\n',mean_grad)

if mean_grad > -1.8
    mean_grad = -1.8;

end

figure;
semilogx(gf_freqs,gf_mag)
xlabel('Frequency (Hz)');
ylabel('Gradient (dB per octave)');
title('Local gradient of selected region');
hold on
semilogx(freq,mag_varSmooth)
grid on;


%% --------------------- Flat Region Design with Two Slopes ---------------

flat_start = 200;
flat_end   = 1000;

fprintf('flat start frequency is %.1f Hz\n', flat_start)

idx_flat_start = find(freq >= flat_start, 1, 'first');
idx_flat_end   = find(freq <= flat_end, 1, 'last');

% Set flat region to average
flatRegion = mean(avgSmooth(idx_flat_start:idx_flat_end));

% Initialize with flat line
L_int = flatRegion * ones(size(freq));

% =========================================================================
% TWO-SLOPE DESIGN (building from flat_start downward)
% =========================================================================
transition_freq = 90;  % Transition point between the two slopes
LF_grad = -3.5;
MF_grad = -8;
% Convert gradients to dB per decade
LF_m_per_decade = LF_grad * log2(10);
MF_m_per_decade = MF_grad * log2(10);

% Apply MF slope (80-200 Hz) - working downward from flat_start
idx_MF = find(freq >= transition_freq & freq < flat_start);
L_int(idx_MF) = flatRegion + MF_m_per_decade .* (log10(freq(idx_MF)) - log10(flat_start));

% Calculate the level at transition point for continuity
level_at_transition = flatRegion + MF_m_per_decade * (log10(transition_freq) - log10(flat_start));

% Apply LF slope (0-80 Hz) - working downward from transition point
idx_LF = find(freq < transition_freq);
L_int(idx_LF) = level_at_transition + LF_m_per_decade .* (log10(freq(idx_LF)) - log10(transition_freq));

L_combined = L_int;
Avg_sum = L_combined;

figure
semilogx(freq, Avg_sum, 'LineWidth', 2, 'DisplayName', 'Target Curve')
grid on
hold on
semilogx(freq, mag, 'DisplayName', 'Measured Response')
xline(flat_start, 'r--', 'DisplayName', 'Flat Start (200 Hz)');
xline(transition_freq, 'g--', 'DisplayName', 'Transition (80 Hz)');
xline(f_cross(1), 'b--', 'DisplayName', 'F cross');
legend('Location', 'best')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Two-Slope Target Curve Design')
xlim([20 500])

%%  LF  roll-off

% Choose magnitude curve for fitting
mag_fit_all = mag_varSmooth;   % or avgSmooth
hi_ratios = 0.7:0.05:1.0;      % f_fit_hi = Fc * ratio
lo_octaves = 1:0.2:2.0;        % span below hi in octaves
best_err = Inf;
best_p = [];
best_range = [];
for r_hi = hi_ratios
    f_fit_hi = f_zeroGrad * r_hi;
    for span_oct = lo_octaves
        f_fit_lo = f_fit_hi / 2^span_oct;
        f_fit_lo = max(30, f_fit_lo);  % enforce 20 Hz minimum

        % Check minimum frequency span
        if (f_fit_hi - f_fit_lo) < 8  % <-- ADD THIS: minimum 15 Hz span
            continue;
        end

        % Extract region
        idx_fit = freq >= f_fit_lo & freq <= f_fit_hi;
        % Require minimum number of points
        if sum(idx_fit) < 10
            continue;
        end
        f_fit  = freq(idx_fit);
        mag_fit = mag_fit_all(idx_fit);
        % Log-frequency axis
        x = log2(f_fit);
        y = mag_fit;
        % Linear fit
        p = polyfit(x, y, 1);
        % Predicted magnitude
        y_hat = polyval(p, x);
        % RMS error
        err = sqrt(mean((y - y_hat).^2));
        % Keep best
        if err < best_err
            best_err   = err;
            best_p     = p;
            best_range = [f_fit_lo, f_fit_hi];
        end
    end
end
% Extract best result
slope_dB_per_oct = best_p(1);
intercept_dB     = best_p(2);
f_fit_lo = best_range(1);
f_fit_hi = best_range(2);
fprintf('Best LF roll-off: %.2f dB/oct\n', slope_dB_per_oct);
fprintf('Fit range: %.1f Hz – %.1f Hz\n', f_fit_lo, f_fit_hi);
fprintf('RMS error: %.2f dB\n', best_err);
figure('Name','Optimal LF Roll-off Fit','NumberTitle','off');
semilogx(freq, mag_fit_all, 'k:', 'DisplayName','Smoothed response');
hold on;
idx_fit = freq >= f_fit_lo & freq <= f_fit_hi;
f_fit = freq(idx_fit);
semilogx(f_fit, polyval(best_p, log2(f_fit)), 'r', 'LineWidth',1.5, ...
    'DisplayName', sprintf('Fit: %.1f dB/oct', slope_dB_per_oct));
xline(f_zeroGrad, '--', 'Fc');
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend show;
xlim([20 f_zeroGrad*1.2]);
% Plot the -minus dB crossing points (frequencies where smoothed == target_dB)
if exist('f_cross','var') && ~isempty(f_cross)

    % highlight the first crossing
    semilogx(f_cross(1), target_dB, 'ro', 'MarkerSize',4, ...
        'DisplayName', sprintf('Selected -%.fdB at %.1f Hz', minus, f_cross(1)));
end
%% blend

% Save the original target before blending
L_combined_original = L_combined;

% Build extended LF slope line from measured roll-off fit
L_LF_slope = polyval(best_p, log2(freq));

% Find intersection
diff_curve = L_combined_original - L_LF_slope;
cross_idx = find(diff(sign(diff_curve)));

if isempty(cross_idx)
    warning('No intersection found between LF slope and target curve.');
    f_intersect = NaN;
else
    i = cross_idx(1);
    f1 = freq(i);
    f2 = freq(i+1);
    d1 = diff_curve(i);
    d2 = diff_curve(i+1);
    f_intersect = f1 + (f2 - f1) * (0 - d1) / (d2 - d1);
end

fprintf('LF slope intersects target at %.1f Hz\n', f_intersect);

% ADJUSTABLE BLEND WITH BIAS CONTROL

% ========== BLEND CONTROL PARAMETERS ==========
blend_width_oct = 1/12;   % Blend width in octaves (1/6 = narrow, 1/2 = wide)
blend_bias = 1.3;        % 0 = follow LF slope more, 1 = follow target more, 0.5 = equal
% ===============================================

% Save the original target before blending
L_combined_original = L_combined;

% Build extended LF slope line from measured roll-off fit
L_LF_slope = polyval(best_p, log2(freq));

% Find intersection
diff_curve = L_combined_original - L_LF_slope;
cross_idx = find(diff(sign(diff_curve)));

if isempty(cross_idx)
    warning('No intersection found between LF slope and target curve.');
    f_intersect = NaN;
else
    i = cross_idx(1);
    f1 = freq(i);
    f2 = freq(i+1);
    d1 = diff_curve(i);
    d2 = diff_curve(i+1);
    f_intersect = f1 + (f2 - f1) * (0 - d1) / (d2 - d1);
end

fprintf('LF slope intersects target at %.1f Hz\n', f_intersect);

L_final = L_combined_original;

if ~isnan(f_intersect)
    % Define blend region
    f_blend_lo = f_intersect * 2^(-blend_width_oct);
    f_blend_hi = f_intersect * 2^( blend_width_oct);

    idx_LF    = freq < f_blend_lo;
    idx_blend = freq >= f_blend_lo & freq <= f_blend_hi;
    idx_HF    = freq > f_blend_hi;

    % Below blend → pure measured LF slope
    L_final(idx_LF) = L_LF_slope(idx_LF);

    % Above blend → keep original target
    L_final(idx_HF) = L_combined_original(idx_HF);

    % ADJUSTABLE BLEND with bias control
    % Normalized position in log space (0 to 1)
    log_pos = (log2(freq(idx_blend)) - log2(f_blend_lo)) / ...
        (log2(f_blend_hi) - log2(f_blend_lo));

    % Apply bias to shift the crossfade point
    % bias < 0.5: favors LF slope longer
    % bias = 0.5: symmetric blend
    % bias > 0.5: transitions to target sooner
    log_pos_biased = log_pos.^(1/blend_bias);
    log_pos_biased = log_pos_biased / max(log_pos_biased);  % renormalize to [0,1]

    % Smooth cosine transition
    alpha = 0.5 * (1 - cos(pi * log_pos_biased));

    % Blend between LF slope and target
    L_final(idx_blend) = ...
        (1 - alpha) .* L_LF_slope(idx_blend) + ...
        alpha       .* L_combined_original(idx_blend);
end

% PLOT THE BLENDED RESULT
figure('Name', sprintf('Adjustable Blend (bias=%.2f)', blend_bias), 'NumberTitle', 'off');

subplot(2,1,1)
semilogx(freq, mag_varSmooth, 'k:', 'DisplayName','Smoothed response');
hold on;
semilogx(freq, L_combined_original, 'b-', 'DisplayName','Original target (flat + slope)');
semilogx(freq, L_LF_slope, 'r--', 'DisplayName','Measured LF roll-off');
if ~isnan(f_intersect)
    xline(f_intersect, 'g--', sprintf('Intersect: %.1f Hz', f_intersect));
end
title('Before Blending')
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Location', 'best');
xlim([20 500]);
ylim([50 100]);

subplot(2,1,2)
semilogx(freq, mag_varSmooth, 'k:', 'DisplayName','Smoothed response');
hold on;
semilogx(freq, L_LF_slope, 'r--', 'DisplayName','Measured LF roll-off');
semilogx(freq, L_combined_original, 'b--', 'DisplayName','Original target');
semilogx(freq, L_final, 'LineWidth', 2, 'DisplayName',...
    sprintf('Blended (bias=%.2f)', blend_bias));
if ~isnan(f_intersect)
    xline(f_intersect, 'g--', sprintf('Intersect: %.1f Hz', f_intersect));
    xregion(f_blend_lo, f_blend_hi, 'FaceColor', [0.9 0.9 0.9], ...
        'FaceAlpha', 0.3, 'DisplayName', 'Blend region');
end
title(sprintf('Adjustable Blend (width=%.2f oct, bias=%.2f)', blend_width_oct, blend_bias))
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Location', 'best');
xlim([20 500]);
ylim([50 100]);

% Optional: Plot the blend weight function to see the bias effect
if ~isnan(f_intersect) && any(idx_blend)
    figure('Name', 'Blend Weight Function', 'NumberTitle', 'off');
    semilogx(freq(idx_blend), alpha, 'b-', 'LineWidth', 2);
    hold on;
    xline(f_intersect, 'g--', 'Intersect');
    xline(f_blend_lo, 'r--', 'Blend start');
    xline(f_blend_hi, 'r--', 'Blend end');
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Blend Weight (0=LF slope, 1=Target)');
    title(sprintf('Blend Weight with Bias=%.2f', blend_bias));
    ylim([0 1]);
end

% Update L_combined with blended result
L_combined = L_final;

%% Condition Input Data

fineSmooth   = mag_varSmooth;   % Variable fine smoothing (1/48 → 1/3 octave)
coarseSmooth = SM(:,5);         % e.g. 1/6-oct smoothing for dips
target       = L_combined;      % Your target slope line

% Compute the difference between smoothed curve and target
diffMag = fineSmooth - target;

% Create a smooth transition weight between 0 and 1
% w = 1 means use fine smoothing, w = 0 means use coarse smoothing
% The tanh() function gives a gentle transition
transitionWidth = 1.5; % in dB, adjust for smoother or sharper transition
w = 0.5 * (1 + tanh(diffMag / transitionWidth));

% Blend between fine and coarse smoothing
mixedSmooth = w .* fineSmooth + (1 - w) .* coarseSmooth;


baseline_smooth = mag;
ref_curve = SM(:,5);   % smoothed reference

% Parameters
min_dip_depth   = 0.5;          % minimum depth to process (dB)
search_range    = [30 300];     % frequency range to analyze dips
smoothing_factor = 1;         % octave smoothing within dip bandwidth
bandwidth_scale = 2;


% Initialize output
mag_dip_conditioned = baseline_smooth;

% Restrict search region
idx_search = (freq >= search_range(1) & freq <= search_range(2));
freq_sub   = baseline_smooth(idx_search);
freq_vec   = freq(idx_search);

% Find local minima (dips) and maxima (peaks)
local_min_idx = [];
local_max_idx = [];

for i = 2:length(freq_sub)-1

    % --- Local minimum (dip)
    if freq_sub(i) < freq_sub(i-1) && freq_sub(i) < freq_sub(i+1)
        left_max  = max(freq_sub(max(1,i-10):i-1));
        right_max = max(freq_sub(i+1:min(end,i+10)));
        depth = min(left_max, right_max) - freq_sub(i);

        if depth > min_dip_depth
            local_min_idx(end+1) = i;
        end
    end

    % --- Local maximum (peak)
    if freq_sub(i) > freq_sub(i-1) && freq_sub(i) > freq_sub(i+1)
        left_min  = min(freq_sub(max(1,i-10):i-1));
        right_min = min(freq_sub(i+1:min(end,i+10)));
        height = freq_sub(i) - max(left_min, right_min);

        if height > min_dip_depth
            local_max_idx(end+1) = i;
        end
    end
end

% Dip mask
dip_mask = zeros(size(freq));

% === PROCESS DIPS
dip_info  = struct([]);
dip_count = 0;

for k = 1:length(local_min_idx)

    idx_local = local_min_idx(k);
    f_center  = freq_vec(idx_local);
    mag_center = freq_sub(idx_local);

    % Index in full-resolution arrays
    idx_center = find(freq == f_center, 1);

    if isempty(idx_center)
        continue;
    end

    ref_center = ref_curve(idx_center);
    dip_depth  = ref_center - mag_center;

    if dip_depth < min_dip_depth
        continue;
    end

    % --- Find left intersection (mag crosses ref_curve)
    iL = idx_center;
    while iL > 1 && mag(iL) < ref_curve(iL)
        iL = iL - 1;
    end

    % --- Find right intersection
    iR = idx_center;
    while iR < length(freq) && mag(iR) < ref_curve(iR)
        iR = iR + 1;
    end

    % Safety
    if iL == idx_center || iR == idx_center || iL >= iR
        continue;
    end

    fL = freq(iL);
    fR = freq(iR);
    BW = fR - fL;

    % --- Distance to nearest peak
    min_dist_to_peak = inf;
    for p = 1:length(local_max_idx)
        min_dist_to_peak = min(min_dist_to_peak, ...
            abs(freq_vec(local_max_idx(p)) - f_center));
    end

    % % --- Q and smoothing width
    % BW = fR - fL;
    % Q  = min(f_center / max(BW,1), 15);
    Q  = f_center / max(BW,1);
    % smooth_width_oct = smoothing_factor * (Q/2);
    % Q_eff = min(Q, 6);
    % smooth_width_oct = bandwidth_scale * smoothing_factor * (Q / 2);

    % % Cap scales with bandwidth control
    % smooth_width_oct = max(0.1, min(smooth_width_oct, 0.5 * bandwidth_scale));
    BW_oct = log2(fR / fL);    % true octave bandwidth
    smooth_width_oct = bandwidth_scale * smoothing_factor * BW_oct;

    smooth_width_oct = max(0.15, min(smooth_width_oct, 1.2));

    if min_dist_to_peak < 2*BW
        smooth_width_oct = smooth_width_oct * 0.5;
    end

    % === STORE DIP (NO HOLES)
    dip_count = dip_count + 1;

    dip_info(dip_count).f_center = f_center;
    dip_info(dip_count).depth = dip_depth;
    dip_info(dip_count).Q = Q;
    dip_info(dip_count).BW = BW;
    dip_info(dip_count).fL = fL;
    dip_info(dip_count).fR = fR;
    % dip_info(dip_count).shoulder_level = shoulder_level;
    dip_info(dip_count).smooth_width = smooth_width_oct;

    fprintf('  Dip %d: f=%.1f Hz, depth=%.2f dB, BW=%.1f Hz, Q=%.2f\n', ...
        dip_count, f_center, dip_depth, BW, Q);

    % Mask and smoothing
    idx_dip = (freq >= fL & freq <= fR);
    dip_mask(idx_dip) = 1;
    %
    % f_smooth_lo = fL * 0.75;
    % f_smooth_hi = fR * 1.25;
    region_scale = 1 + 0.2*(bandwidth_scale - 1);  % gentle expansion
    f_smooth_lo = fL / region_scale;
    f_smooth_hi = fR * region_scale;

    idx_smooth  = (freq >= f_smooth_lo & freq <= f_smooth_hi);

    freq_s = freq(idx_smooth);
    mag_s  = mag(idx_smooth);

    width = smooth_width_oct;
    smoothed_vals = zeros(size(freq_s));

    for i = 1:length(freq_s)
        dlog = log2(freq_s / freq_s(i));
        w = exp(-(dlog.^2)/(2*width^2));
        w = w / sum(w);
        smoothed_vals(i) = sum(w .* mag_s);
    end

    % blend_margin = 0.05;
    blend_margin = 0.05 * bandwidth_scale;
    blend_margin = min(blend_margin, 0.2);  % safety cap

    f_blend_lo = fL*(1-blend_margin);
    f_blend_hi = fR*(1+blend_margin);

    for i = 1:length(freq)
        if freq(i) >= f_blend_lo && freq(i) <= f_blend_hi
            if freq(i) >= fL && freq(i) <= fR
                bw = 1;
            elseif freq(i) < fL
                bw = (freq(i)-f_blend_lo)/(fL-f_blend_lo);
            else
                bw = (f_blend_hi-freq(i))/(f_blend_hi-fR);
            end
            bw = max(0,min(1,0.5*(1-cos(pi*bw))));

            [~,si] = min(abs(freq_s - freq(i)));
            if bw > 0 && smoothed_vals(si) > baseline_smooth(i)
                mag_dip_conditioned(i) = ...
                    (1-bw)*mag_dip_conditioned(i) + bw*smoothed_vals(si);
            end
        end
    end


    % % conformation plots
    % figure;
    % semilogx(freq, mag, 'k'); hold on;
    % semilogx(freq, ref_curve, 'b', 'LineWidth', 1.5);
    % xline(fL,'r'); xline(f_center,'r--'); xline(fR,'r');
    % grid on;
    % title(sprintf('Dip @ %.1f Hz, depth %.2f dB', f_center, dip_depth));


end

fprintf('Processed %d valid dips\n', dip_count);

mag_final_conditioned = mag_dip_conditioned;


figure('Name', 'Dip-ONLY Smoothing with Peak Protection', 'NumberTitle', 'off');

% Plot 1: Show dip and peak detection
subplot(3,1,1);
semilogx(freq, mag, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Original data');
hold on;

semilogx(freq, SM(:,5),'LineWidth', 0.5, ...
    'DisplayName', 'Variable smoothing (baseline)');
semilogx(freq, mag_varSmooth, 'b-', 'DisplayName', 'Variable smoothing');
semilogx(freq, mag_final_conditioned, 'r-', 'LineWidth', 1.5, ...
    'DisplayName', '1/3rd oct smoothing');
% Mark detected dips (red)
for k = 1:length(dip_info)
    xregion(dip_info(k).fL, dip_info(k).fR, ...
        'FaceColor', [1 0.8 0.8], 'FaceAlpha', 0.5, ...
        'HandleVisibility', 'off');
    % plot(dip_info(k).f_center, ...
    %      dip_info(k).shoulder_level - dip_info(k).depth, ...
    %      'rv', 'MarkerSize', 8, 'MarkerFaceColor', 'r', ...
    %      'DisplayName', 'Dips (will smooth)');
end

% Mark detected peaks (green)
for k = 1:length(local_max_idx)
    f_peak = freq_vec(local_max_idx(k));
    mag_peak = freq_sub(local_max_idx(k));
    plot(f_peak, mag_peak, '^', 'MarkerSize', 8, ...
        'Color', [0 0.7 0], 'MarkerFaceColor', [0 0.7 0], ...
        'DisplayName', 'Peaks (protected)');
end

% Show search range
xline(search_range(1), 'k--', '30 Hz', 'HandleVisibility', 'off');
xline(search_range(2), 'k--', '300 Hz', 'HandleVisibility', 'off');

xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Dip and Peak Detection (30-300 Hz only)');
legend('Location', 'best');
grid on;
xlim([20 500]);

% Plot 2: Show smoothing effect
subplot(3,1,2);
semilogx(freq, mag_varSmooth, 'b-', 'LineWidth', 0.5, ...
    'DisplayName', 'Variable smoothing (baseline)');
hold on;
semilogx(freq, SM(:,5),'LineWidth', 0.5, ...
    'DisplayName', 'Variable smoothing (baseline)');
semilogx(freq, mag_final_conditioned, 'r-', 'LineWidth', 1.5, ...
    'DisplayName', '1/3rd oct smoothing');

% Shade dip regions that were smoothed
for k = 1:length(dip_info)
    xregion(dip_info(k).fL, dip_info(k).fR, ...
        'FaceColor', [1 0.8 0.8], 'FaceAlpha', 0.3, ...
        'HandleVisibility', 'off');
end

% Shade protected peak regions
for k = 1:length(local_max_idx)
    f_peak = freq_vec(local_max_idx(k));
    protect_width = 1/6;
    f_protect_lo = f_peak * 2^(-protect_width);
    f_protect_hi = f_peak * 2^(protect_width);
    xregion(f_protect_lo, f_protect_hi, ...
        'FaceColor', [0.8 1 0.8], 'FaceAlpha', 0.3, ...
        'HandleVisibility', 'off');
end

xline(search_range(1), 'k--', 'HandleVisibility', 'off');
xline(search_range(2), 'k--', 'HandleVisibility', 'off');

xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Final Result: Dips Smoothed, Peaks Preserved');
legend('Location', 'best');
grid on;
xlim([20 500]);

% Plot 3: Show the difference
subplot(3,1,3);
difference = mag_final_conditioned - mag_varSmooth;
semilogx(freq, difference, 'k-', 'LineWidth', 1.5, ...
    'DisplayName', 'Change (final - baseline)');
hold on;
yline(0, 'b--', 'HandleVisibility', 'off');

% Shade regions
for k = 1:length(dip_info)
    xregion(dip_info(k).fL, dip_info(k).fR, ...
        'FaceColor', [1 0.8 0.8], 'FaceAlpha', 0.3, ...
        'DisplayName', 'Dip regions');
end

for k = 1:length(local_max_idx)
    f_peak = freq_vec(local_max_idx(k));
    protect_width = 1/6;
    f_protect_lo = f_peak * 2^(-protect_width);
    f_protect_hi = f_peak * 2^(protect_width);
    xregion(f_protect_lo, f_protect_hi, ...
        'FaceColor', [0.8 1 0.8], 'FaceAlpha', 0.3, ...
        'DisplayName', 'Protected peaks');
end

xline(search_range(1), 'k--', 'HandleVisibility', 'off');
xline(search_range(2), 'k--', 'HandleVisibility', 'off');

xlabel('Frequency (Hz)');
ylabel('Change (dB)');
title('Difference: Positive = Dip Filled, Zero at Peaks');
legend('Location', 'best');
grid on;
xlim([20 500]);

% Use conditioned data for room correction
mixedSmooth = mag_final_conditioned;

figure('Name','Mixed Variable Smoothing');
semilogx(freq, mag_varSmooth, 'Color',[0.7 0.7 0.7],'DisplayName','varable smoothing');
hold on;
semilogx(freq, fineSmooth,'DisplayName','Fine smoothing (variable)');
semilogx(freq, coarseSmooth,'DisplayName','Coarse smoothing (1/6 octave)');
semilogx(freq, mag, 'b','LineWidth',1.5,'DisplayName','Input data');
semilogx(freq, mag_final_conditioned, 'k','LineWidth',1.5,'DisplayName','DA smooting');
semilogx(freq, target, 'g--','DisplayName','Target line');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Mixed Smoothing: Preserve Peaks, Smooth Dips');
legend show; grid on;
xlim([30 1000]); ylim([30 max(mag)+10])


%% Peaking filters with Frequency-Dependent Limiting - Two Phase Approach
%
% Avg_sum = L_combined;
% maxIters = 15;
% peak_store1 = [];
% peakID = 0;
%
% % correction_range = [f_cross(1,1) 300];
% correction_range = [30 300];
% sos_all = [];
% mag_current = mixedSmooth;
%
% limitConfig = generateLimitBands(30, 500, 8, ...
%     'boost_start', 20, 'boost_end', 20,'cut_start',20,'cut_end',20);
%
% % PHASE 1: FILL DIPS
% fprintf('\n=== PHASE 1: Filling Dips ===\n');
%
% for iter = 1:maxIters
%     % Calculate residual
%     resid = mag_current - Avg_sum(:);
%
%     % === FIND WORST DIP (inline) ===
%     thresh = -1.0;  % minimum dip depth (dB)
%
%     % Restrict search range
%     idx_range = (freq >= correction_range(1) & freq <= correction_range(2));
%     resid_sub = resid(idx_range);
%     freq_sub = freq(idx_range);
%
%     isDip = resid_sub < thresh;
%     d = diff([0; isDip; 0]);
%     startIdx = find(d == 1);
%     endIdx = find(d == -1) - 1;
%
%     if isempty(startIdx)
%         fprintf('No more significant dips found after %d iterations.\n', iter-1);
%         break;
%     end
%
%     % Analyze each dip region
%     nDips = numel(startIdx);
%     dipMeans = zeros(nDips,1);
%     dipCenters = zeros(nDips,1);
%     dipWidths = zeros(nDips,1);
%     dipImportance = zeros(nDips,1);
%
%     for k = 1:nDips
%         idx = startIdx(k):endIdx(k);
%         dipMeans(k) = mean(resid_sub(idx));        % average depth (dB)
%         dipCenters(k) = mean(freq_sub(idx));       % center freq (Hz)
%         dipWidths(k) = freq_sub(endIdx(k)) - freq_sub(startIdx(k)); % bandwidth (Hz)
%
%         dipDepth = abs(min(resid_sub(idx)));
%         f1 = freq_sub(startIdx(k));
%         f2 = freq_sub(endIdx(k));
%         width_oct = log2(f2 / f1);
%         dipImportance(k) = width_oct / dipDepth;  % importance metric
%     end
%
%     % Pick most significant dip
%     [~, idxBest] = max(dipImportance);
%     f0 = dipCenters(idxBest);
%     amp = dipMeans(idxBest);
%     width = dipWidths(idxBest);
%
%     % Calculate Q for dip
%     Q_base = f0 / max(width, 1e-6);
%     raw_Q = max(2, min(10, Q_base));  % clamp Q to safe range
%
%     % Extract dip parameters
%     raw_gain = abs(amp);  % positive gain to boost dip
%
%
%     % === END FIND WORST DIP ===
%
%     % Apply frequency limits
%     [gain, Qval, band_idx] = applyFrequencyLimits(f0, raw_gain, raw_Q, limitConfig);
%
%     % Design and apply filter
%     sos = designPeakingBiquad(f0, fs, gain, Qval);
%     mag_current = applyFilter(mag_current, freq, fs, sos);
%
%     % Store filter data
%     peakID = peakID + 1;
%     peak_store1(end+1,:) = [f0, gain, Qval, peakID];
%     sos_all = [sos_all; sos];
%
%     % Plot progress
%     plotIterationProgress(iter, freq, mixedSmooth, mag_current, Avg_sum, sos_all);
% end
%

% === PARAMETRIC EQ CORRECTION STAGE WITH DIP PROTECTION ===
Avg_sum = target;  % Use the same target
maxIters = 15;
peak_store1 = [];
peakID = 0;

correction_range = [30 300];
sos_all = [];
mag_current = mag_final_conditioned;  % Start from conditioned data

limitConfig = generateLimitBands(30, 500, 8, ...
    'boost_start', 20, 'boost_end', 20,'cut_start',20,'cut_end',20);

% === CREATE DIP PROTECTION MAP ===
% Build a gain limit array based on detected dips
max_boost_limit = ones(size(freq)) * 20;  % default: allow up to 20 dB boost

for k = 1:length(dip_info)
    % Get dip parameters
    f_center = dip_info(k).f_center;
    fL = dip_info(k).fL;
    fR = dip_info(k).fR;
    dip_depth = dip_info(k).depth;

    % Calculate protection zone (wider than the dip itself)
    protection_scale = 1.5;  % extend protection 50% beyond dip edges
    f_protect_lo = fL / protection_scale;
    f_protect_hi = fR * protection_scale;

    % Find indices in protection zone
    idx_protect = (freq >= f_protect_lo & freq <= f_protect_hi);

    % Limit boost based on dip depth
    % Deep dips (>6 dB) = very limited boost (2-3 dB max)
    % Moderate dips (3-6 dB) = moderate boost (4-6 dB max)
    % Shallow dips (<3 dB) = normal boost allowed
    if dip_depth > 6
        allowed_boost = 2.0;  % very limited for deep nulls
    elseif dip_depth > 3
        allowed_boost = 4.0;  % moderate for moderate dips
    else
        allowed_boost = 8.0;  % more freedom for shallow dips
    end

    % Apply taper from center to edges of protection zone
    for i = 1:length(freq)
        if idx_protect(i)
            % Distance from dip center (in octaves)
            dist_oct = abs(log2(freq(i) / f_center));
            BW_oct = log2(fR / fL);

            % Taper: full restriction at center, relaxes toward edges
            if dist_oct < BW_oct/2
                % Inside the dip bandwidth
                taper = 0;  % full restriction
            else
                % In the protection zone beyond dip
                taper = (dist_oct - BW_oct/2) / (log2(protection_scale));
                taper = min(1, max(0, taper));  % clamp 0-1
            end

            % Blend between restricted and unrestricted
            limit_here = allowed_boost + taper * (20 - allowed_boost);

            % Take the minimum (most restrictive)
            max_boost_limit(i) = min(max_boost_limit(i), limit_here);
        end
    end

    fprintf('  Dip protection zone %d: %.1f-%.1f Hz, max boost: %.1f dB\n', ...
        k, f_protect_lo, f_protect_hi, allowed_boost);
end

% === INITIALIZE EXCLUSION MASK ===
% Track regions where we've already tried corrections
exclusion_mask = zeros(size(freq));  % 0 = allowed, 1 = excluded
exclusion_history = [];  % Store [f_center, fL, fR, gain_attempted]

% SINGLE OPTIMIZED PHASE: BOTH CUTS AND BOOSTS WITH DIP PROTECTION
fprintf('\n=== OPTIMIZED CORRECTION PHASE (with dip protection) ===\n');

for iter = 1:maxIters
    % Calculate residual using consistent target variable
    resid = mag_current - target;

    % Define correction range
    idx_range = (freq >= correction_range(1) & freq <= correction_range(2));
    freq_sub = freq(idx_range);
    resid_sub = resid(idx_range);

    % Apply exclusion mask to residual
    exclusion_sub = exclusion_mask(idx_range);
    resid_sub_masked = resid_sub;
    resid_sub_masked(exclusion_sub > 0) = 0;  % Zero out excluded regions

    % Find worst deviation (peak OR dip) in non-excluded regions
    [abs_resid, max_idx] = max(abs(resid_sub_masked));

    if abs_resid < 1.0  % threshold in dB
        fprintf('Residual below threshold after %d iterations.\n', iter-1);
        break;
    end

    % Get deviation parameters
    amp = resid_sub(max_idx);  % positive for peak, negative for dip
    f0 = freq_sub(max_idx);
    is_boost = (amp < 0);  % negative residual means we need to boost

    fprintf('\nIter %d Analysis:\n', iter);
    if is_boost
        fprintf('  Dip at %.1f Hz, needs %.2f dB boost\n', f0, abs(amp));
    else
        fprintf('  Peak at %.1f Hz, needs %.2f dB cut\n', f0, abs(amp));
    end

    % === MEASURE THE ACTUAL DEVIATION WIDTH ===
    [~, peak_idx] = min(abs(freq - f0));

    % Find the bandwidth where deviation actually exists
    deviation_threshold = abs(amp) * 0.1;  % 10% of deviation amplitude

    % Search left from deviation
    iL = max_idx;
    while iL > 1 && abs(resid_sub(iL)) > deviation_threshold
        iL = iL - 1;
    end
    fL_deviation = freq_sub(iL);

    % Search right from deviation
    iR = max_idx;
    while iR < length(freq_sub) && abs(resid_sub(iR)) > deviation_threshold
        iR = iR + 1;
    end
    fR_deviation = freq_sub(iR);

    % Calculate the actual deviation bandwidth
    BW_deviation = fR_deviation - fL_deviation;
    BW_oct_deviation = log2(fR_deviation / fL_deviation);

    % Calculate Q needed to match this bandwidth
    Q_needed = f0 / BW_deviation;
    Q_needed = max(Q_needed, 2);   % minimum Q of 2
    Q_needed = min(Q_needed, 15);  % maximum Q of 15 for stability

    fprintf('  Deviation bandwidth: %.1f-%.1f Hz (%.3f octaves)\n', ...
        fL_deviation, fR_deviation, BW_oct_deviation);
    fprintf('  Required Q: %.2f\n', Q_needed);

    % === CHECK BOOST LIMIT AT THIS FREQUENCY (if boosting) ===
    [~, f0_idx] = min(abs(freq - f0));
    boost_limit_here = max_boost_limit(f0_idx);

    limited_boost = false;
    if is_boost
        required_boost = abs(amp);
        if required_boost > boost_limit_here
            fprintf('  >> Limited from %.2f dB to %.2f dB boost\n', ...
                required_boost, boost_limit_here);
            amp = -boost_limit_here;  % limit the boost
            limited_boost = true;
        end
    end

    % === DESIGN SURGICAL FILTER ===
    gain_needed = -amp * 1.05;  % 5% extra to ensure we clear the target

    % Design the filter
    sos = designPeakingBiquad(f0, fs, gain_needed, Q_needed);

    % === VERIFY FILTER RESPONSE ===
    % Calculate what this filter actually does across frequency
    filter_response = applyFilter(zeros(size(freq)), freq, fs, sos);

    % Find -3dB bandwidth of the filter
    if is_boost
        peak_filter_response = max(filter_response);  % most positive for boost
        half_gain = peak_filter_response / 2;

        iL_filter = peak_idx;
        while iL_filter > 1 && filter_response(iL_filter) > half_gain
            iL_filter = iL_filter - 1;
        end

        iR_filter = peak_idx;
        while iR_filter < length(freq) && filter_response(iR_filter) > half_gain
            iR_filter = iR_filter + 1;
        end
    else
        peak_filter_response = min(filter_response);  % most negative for cut
        half_gain = peak_filter_response / 2;

        iL_filter = peak_idx;
        while iL_filter > 1 && filter_response(iL_filter) < half_gain
            iL_filter = iL_filter - 1;
        end

        iR_filter = peak_idx;
        while iR_filter < length(freq) && filter_response(iR_filter) < half_gain
            iR_filter = iR_filter + 1;
        end
    end

    BW_filter = freq(iR_filter) - freq(iL_filter);
    BW_oct_filter = log2(freq(iR_filter) / freq(iL_filter));

    fprintf('  Filter -3dB bandwidth: %.1f-%.1f Hz (%.3f octaves)\n', ...
        freq(iL_filter), freq(iR_filter), BW_oct_filter);

    % === CHECK COLLATERAL DAMAGE BEFORE APPLYING ===
    mag_test = applyFilter(mag_current, freq, fs, sos);

    % Define "safe zone" - slightly wider than deviation zone
    safe_margin = 1.5;
    f_safe_lo = fL_deviation / safe_margin;
    f_safe_hi = fR_deviation * safe_margin;

    % Check impact outside safe zone
    outside_safe = (freq < f_safe_lo) | (freq > f_safe_hi);
    outside_safe = outside_safe & (freq >= correction_range(1)) & (freq <= correction_range(2));

    collateral_change = abs(mag_test(outside_safe) - mag_current(outside_safe));
    max_collateral = max(collateral_change);
    mean_collateral = mean(collateral_change);

    fprintf('  Collateral damage outside %.1f-%.1f Hz:\n', f_safe_lo, f_safe_hi);
    fprintf('    Max: %.3f dB, Mean: %.3f dB\n', max_collateral, mean_collateral);

    % === CHECK IMPROVEMENT ===
    % Calculate RMS improvement in the deviation region
    idx_deviation = (freq >= fL_deviation & freq <= fR_deviation);
    resid_before = mag_current(idx_deviation) - target(idx_deviation);
    resid_after = mag_test(idx_deviation) - target(idx_deviation);
    rms_before = sqrt(mean(resid_before.^2));
    rms_after = sqrt(mean(resid_after.^2));
    improvement = rms_before - rms_after;

    fprintf('  RMS improvement in target region: %.3f dB\n', improvement);

    % === DECISION: APPLY OR EXCLUDE ===
    should_exclude = false;
    exclude_reason = '';

    % Check 1: Minimal improvement
    if improvement < 0.1
        should_exclude = true;
        exclude_reason = sprintf('minimal improvement (%.3f dB)', improvement);
    end

    % Check 2: For boost filters, if gain is much smaller than needed
    if is_boost && limited_boost && abs(gain_needed) < abs(amp) * 0.3
        should_exclude = true;
        exclude_reason = sprintf('ineffective boost (%.1f dB vs %.1f dB needed)', ...
            abs(gain_needed), abs(amp));
    end

    % Check 3: Excessive collateral damage
    if max_collateral > abs(amp) * 0.4
        fprintf('  WARNING: High collateral damage (%.3f dB vs %.3f dB correction)\n', ...
            max_collateral, abs(amp));
    end

    if should_exclude
        fprintf('  >> Excluding region: %s\n', exclude_reason);

        % Add to exclusion mask
        exclusion_margin = 1.3;
        f_exclude_lo = fL_deviation / exclusion_margin;
        f_exclude_hi = fR_deviation * exclusion_margin;
        idx_exclude = (freq >= f_exclude_lo & freq <= f_exclude_hi);
        exclusion_mask(idx_exclude) = 1;

        % Record in history
        exclusion_history(end+1, :) = [f0, f_exclude_lo, f_exclude_hi, gain_needed];

        fprintf('  >> Skipping filter application, moving to next deviation\n\n');
        continue;
    end
    % === APPLY FILTER ===
    mag_before = mag_current;
    mag_current = mag_test;

    % Verify actual improvement at peak (find the peak in the full freq array)
    [~, peak_idx_full] = min(abs(freq - f0));
    resid_before_full = mag_before - target;
    resid_after_full = mag_current - target;
    actual_improvement_at_peak = abs(resid_before_full(peak_idx_full)) - abs(resid_after_full(peak_idx_full));

    fprintf('  Applied filter: improvement at peak = %.3f dB\n', actual_improvement_at_peak);

    % Store filter data
    peakID = peakID + 1;
    peak_store1(end+1,:) = [f0, gain_needed, Q_needed, peakID];
    sos_all = [sos_all; sos];

    % === DETAILED PLOT WITH CUMULATIVE HISTORY ===
    figure('Name', sprintf('Filter %d - Correction', iter));

    % Main plot - Cumulative effect
    subplot(2,1,1);
    semilogx(freq, target, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Target');
    hold on;

    % Show original conditioned input (iteration 0)
    semilogx(freq, mag_final_conditioned, 'Color', [0.85 0.85 0.85], ...
        'LineWidth', 0.8, 'DisplayName', 'Original (iter 0)');

    % % Show all previous iterations in light gray
    % if iter > 1
    %     % Reconstruct each iteration's magnitude
    %     mag_history = mag_final_conditioned;
    %     for prev_iter = 1:iter-1
    %         sos_prev = sos_all(prev_iter, :);
    %         mag_history = applyFilter(mag_history, freq, fs, sos_prev);
    %
    %         % Plot previous iterations in progressively lighter gray
    %         gray_level = 0.6 + (prev_iter / iter) * 0.2;  % 0.6 to 0.8
    %         semilogx(freq, mag_history, 'Color', [gray_level gray_level gray_level], ...
    %                  'LineWidth', 0.8, 'HandleVisibility', 'off');
    %     end
    % end
    %
    % % Highlight the step from previous to current (this filter's effect)
    % semilogx(freq, mag_before, ':', 'Color', [0.3 0.3 0.9], ...
    %          'LineWidth', 1.5, 'DisplayName', sprintf('Before filter %d', iter));
    %
    % Show current state after this filter (bold)
    semilogx(freq,mag,'DisplayName','input responce')
    semilogx(freq, mag_current, 'LineWidth', 2.2, ...
        'Color', [0 0.45 0.9], 'DisplayName', sprintf('After filter %d', iter));

    % Highlight deviation zone with color based on type
    if is_boost
        xregion(fL_deviation, fR_deviation, 'FaceColor', [1 0.8 0.8], 'FaceAlpha', 0.4, ...
            'DisplayName', 'Dip region');
    else
        xregion(fL_deviation, fR_deviation, 'FaceColor', [0.8 0.8 1], 'FaceAlpha', 0.4, ...
            'DisplayName', 'Peak region');
    end

    % Highlight safe zone
    xregion(f_safe_lo, fL_deviation, 'FaceColor', [1 1 0.8], 'FaceAlpha', 0.2, ...
        'HandleVisibility', 'off');
    xregion(fR_deviation, f_safe_hi, 'FaceColor', [1 1 0.8], 'FaceAlpha', 0.2, ...
        'HandleVisibility', 'off');

    % Show excluded regions (grayed out)
    for h = 1:size(exclusion_history, 1)
        xregion(exclusion_history(h, 2), exclusion_history(h, 3), ...
            'FaceColor', [0.9 0.9 0.9], 'FaceAlpha', 0.6, ...
            'HandleVisibility', 'off');
    end

    % Mark filter location
    xline(f0, 'r--', sprintf('%.0f Hz', f0), 'LineWidth', 1.2, ...
        'LabelVerticalAlignment', 'bottom');

    % Calculate cumulative RMS error
    cumulative_resid = mag_current - target;
    cumulative_rms = sqrt(mean(cumulative_resid(idx_range).^2));

    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title(sprintf('Cumulative: %d filters | This: f0=%.0f Hz, Q=%.1f, Gain=%+.1f dB | RMS: %.2f dB', ...
        iter, f0, Q_needed, gain_needed, cumulative_rms));
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    xlim([20 1000]);

    % Filter response plot
    subplot(2,1,2);
    semilogx(freq,mag_current -mag_final_conditioned, 'LineWidth',2.2)
    xlim([20 1000]);
    grid minor;
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');

    title('total cumulative residual')

    %
    % % Plot all previous filter responses in light colors
    % if iter > 1
    %     for prev_iter = 1:iter-1
    %         sos_prev = sos_all(prev_iter, :);
    %         filter_resp_prev = applyFilter(zeros(size(freq)), freq, fs, sos_prev);
    %
    %         % Color code: red for cuts, blue for boosts
    %         f0_prev = peak_store1(prev_iter, 1);
    %         gain_prev = peak_store1(prev_iter, 2);
    %
    %         if gain_prev > 0
    %             color_prev = [0.8 0.9 1];  % light blue for boost
    %         else
    %             color_prev = [1 0.9 0.9];  % light red for cut
    %         end
    %
    %         semilogx(freq, filter_resp_prev, 'Color', color_prev, ...
    %                  'LineWidth', 0.8, 'HandleVisibility', 'off');
    %         hold on;
    %     end
    % end
    %
    % % Plot current filter response (bold)
    % if is_boost
    %     filter_color = [0 0.5 1];  % blue for boost
    % else
    %     filter_color = [1 0 0];    % red for cut
    % end
    %
    % semilogx(freq, filter_response, 'Color', filter_color, 'LineWidth', 2, ...
    %          'DisplayName', sprintf('Filter %d response', iter));
    % hold on;
    %
    % % Mark half-gain line
    % yline(peak_filter_response/2, 'k--', '-3dB', 'LineWidth', 0.8, ...
    %       'DisplayName', '-3dB line', 'LabelHorizontalAlignment', 'left');
    %
    % % Highlight filter bandwidth
    % if is_boost
    %     xregion(freq(iL_filter), freq(iR_filter), 'FaceColor', [0.8 0.9 1], ...
    %             'FaceAlpha', 0.4, 'DisplayName', 'Filter -3dB BW');
    % else
    %     xregion(freq(iL_filter), freq(iR_filter), 'FaceColor', [1 0.9 0.9], ...
    %             'FaceAlpha', 0.4, 'DisplayName', 'Filter -3dB BW');
    % end
    %
    % % Show deviation zone for comparison
    % if is_boost
    %     xregion(fL_deviation, fR_deviation, 'FaceColor', [1 0.8 0.8], 'FaceAlpha', 0.2, ...
    %             'DisplayName', 'Problem zone');
    % else
    %     xregion(fL_deviation, fR_deviation, 'FaceColor', [0.8 0.8 1], 'FaceAlpha', 0.2, ...
    %             'DisplayName', 'Problem zone');
    % end
    %
    % xline(f0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    %
    % % Add text annotation showing all filters applied so far
    % filter_summary = sprintf('Filters 1-%d shown', iter);
    % text(0.02, 0.98, filter_summary, 'Units', 'normalized', ...
    %      'VerticalAlignment', 'top', 'FontSize', 9, ...
    %      'BackgroundColor', [1 1 1 0.8]);
    %
    % xlabel('Frequency (Hz)');
    % ylabel('Filter Gain (dB)');
    % title(sprintf('Individual Filter: BW=%.1f Hz (%.3f oct), Q=%.1f | Total: %d filters', ...
    %               BW_filter, BW_oct_filter, Q_needed, iter));
    % legend('Location', 'best', 'FontSize', 8);
    % grid on;
    % xlim([20 1000]);
    %
    % % Add y-axis limits based on filter type
    % if is_boost
    %     ylim([min(filter_response)-1, max(filter_response)+1]);
    % else
    %     ylim([min(filter_response)-1, max(filter_response)+1]);
    % end
    %
    drawnow;
    pause(0.3);  % Brief pause to see the progression

    fprintf('Iter %d: f0=%.1f Hz, gain=%.2f dB, Q=%.2f, improvement=%.3f dB, cumulative RMS=%.3f dB\n', ...
        iter, f0, gain_needed, Q_needed, improvement, cumulative_rms);
end

fprintf('\nTotal filters applied: %d\n', size(sos_all, 1));
fprintf('Excluded regions: %d\n', size(exclusion_history, 1));

PPB = mag_current;
num_dip_filters = size(sos_all, 1);
PPB_dips_filled = mag_current;


%
% for iter = 1:15-num_dip_filters
%     resid_p = mag_current - L_combined;
%     thresh_p = 1;%dB
%     idx_range = (freq >= correction_range(1) & freq <= correction_range(2));
%     freq_sub = freq(idx_range);
%     resid_sub = resid_p(idx_range);
%     [pk, lc] = s_findpeaks(resid_sub, thresh, 10);
%     [~, idxMax] = max(pk);
%     amp = pk(idxMax);
%     idx = lc(idxMax);
%     f0 = freq_sub(idx);
%
%     % find Q
%     half_height = amp *0.7;
%     % Search left from peak
%     iL = idx;
%     while iL > 1 && resid_sub(iL) > half_height
%         iL = iL - 1;
%     end
%     fL = freq_sub(iL);
%     % Search right from peak
%     iR = idx;
%     while iR < length(freq_sub) && resid_sub(iR) > half_height
%         iR = iR + 1;
%     end
%     fR = freq_sub(iR);
%     % Calculate bandwidth and Q
%     BW = fR - fL;
%     Q_initial = (f0 / BW);
%     Q_initial = min(Q_initial, 10);  % Cap initial Q at 15
%
%     % ------------------------- optimisation ------------------------------
%
%     % Define optimization region (around the peak)
%     opt_range = [fL-15, fR+15];  % optimize within the bandwidth
%     opt_idx = (freq >= opt_range(1) & freq <= opt_range(2));
%
%     % Initial guess
%     f0_init = f0;
%     gain_init = -amp;
%     Q_init = Q_initial;
%
%     % Objective function: minimize RMS residual in the optimization region
%     objective = @(params) optimizeFilter(params, mag_current, L_combined, freq, fs, opt_idx);
%
%     % Parameter bounds: [f0, gain, Q]
%     lb = [f0*0.5,  -amp*2,  Q_initial*0.2];     % lower bounds
%     ub = [f0*1.5,  -amp*0.2,  min(Q_initial*2, 15)];  % upper bound capped at 15
%     x0 = [f0_init, gain_init, Q_init];  % initial guess
%
%     % Run optimization with bounds using fmincon
%     options = optimoptions('fmincon', 'Display', 'off', 'TolX', 0.1, 'TolFun', 0.01);
%     [params_opt, fval] = fmincon(objective, x0, [], [], [], [], lb, ub, [], options);
%
%     % Extract optimized parameters
%     f0_opt = params_opt(1);
%     gain_opt = params_opt(2);
%     Q_opt = params_opt(3);
%     Q_opt = min(Q_opt, 15);  % Extra safety check to ensure Q ≤ 15
%
%     % generate optimized filter
%     sos = designPeakingBiquad(f0_opt, fs, gain_opt, Q_opt);
%
%     % Apply the filter
%     mag_after = applyFilter(mag_current, freq, fs, sos);
%
%     % Calculate new residual to check
%     resid_after = mag_after - L_combined;
%
%     % === PLOT THE RESULTS ===
%     figure('Name', 'Filter Result - Optimized');
%
%     % Top plot: magnitude response
%     % subplot(2,1,1);
%     semilogx(freq, mag_current, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Before filter');
%     hold on;
%     semilogx(freq, mag_after, 'r-', 'LineWidth', 1.5, 'DisplayName', 'After filter');
%     semilogx(freq, L_combined, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Target');
%     semilogx(freq, L_combined+3, '--', 'DisplayName', 'Target + 3','Color',[0.5 0.5 0.5]);
%     semilogx(freq, L_combined-3, '--', 'DisplayName', 'Target - 3','Color',[0.5 0.5 0.5]);
%
%     xline(f0_opt, 'k--', sprintf('%.1f Hz', f0_opt));
%     xlabel('Frequency (Hz)');
%     ylabel('Magnitude (dB)');
%     title(sprintf('Optimized Filter: f0=%.1f Hz, gain=%.2f dB, Q=%.2f', f0_opt, gain_opt, Q_opt));
%     legend show;
%     grid on;
%     xlim([20 1000]);
%
%     mag_current = mag_after;
%
%     peakID = peakID + 1;
%     peak_store1(end+1,:) = [f0_opt, gain_opt, Q_opt, peakID];
%     sos_all = [sos_all; sos];
%
% end
%
% PPB = mag_current;
%
% % optimisation helper function
% function rms_error = optimizeFilter(params, mag_current, target, freq, fs, opt_idx)
% % Extract parameters
% f0 = params(1);
% gain = params(2);
% Q = params(3);
%
% % Design filter with these parameters
% sos = designPeakingBiquad(f0, fs, gain, Q);
%
% % Apply filter
% mag_after = applyFilter(mag_current, freq, fs, sos);
%
% % Calculate residual in optimization region only
% resid = mag_after - target;
% resid_opt = resid(opt_idx);
%
% % RMS error
% rms_error = sqrt(mean(resid_opt.^2));
% end


%%  expected responce without limiters
mag_expected = mag;
for i = 1:size(sos_all, 1)
    sos_i = sos_all(i, :);
    mag_expected = applyFilter(mag_expected, freq, fs, sos_i);
end
figure('Name','expected responce without limits')
% subplot(2,1,1)
semilogx(freq, mag, 'Color',[0.5 0.5 0.5],'LineWidth', 1, 'DisplayName', 'REW Input');
hold on;
semilogx(freq, mag_expected, 'LineWidth',1.5,'DisplayName', 'Mag with applied filters');
semilogx(freq, L_combined, 'g--','LineWidth',1.5,'DisplayName', 'Target');
% semilogx(freq, L_combined+3, '--','DisplayName', 'Target + 3','Color',[0.8 0.5 0.5]);
% semilogx(freq, L_combined-3, '--','DisplayName', 'Target - 3','Color',[0.8 0.5 0.5]);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title(sprintf('Expected Response - %d Filters Applied', size(sos_all, 1)));
legend show;
grid minor;
xlim([20 1000]);


user_slope_dB_per_oct = [1.2];   % e.g. 3 or [] to auto-compute to hit limit_max at f_cap

H_plot = (mixedSmooth - mag_current)*-1;


% Compute slope:
if isempty(user_slope_dB_per_oct)

    slope_dB_per_oct = limit_max / log2(f_cap / f_ref);
else
    slope_dB_per_oct = user_slope_dB_per_oct;
end

limit_curve = 3*(ones(size(freq)));
idx_above_ref = freq > f_ref;
limit_curve(idx_above_ref) = 3 + slope_dB_per_oct .* log2(freq(idx_above_ref) ./ f_ref);


limit_curve(limit_curve < 0) = 2;
limit_curve(limit_curve > limit_max) = limit_max;

figure('Name','Limiter Rising Envelope');
semilogx(freq, limit_curve, 'LineWidth', 1.6, 'DisplayName', 'Limiter Envelope');
hold on;
semilogx(freq, H_plot, '--', 'DisplayName', 'Applied Boost');
xlabel('Frequency (Hz)');
ylabel('Boost (dB)');
grid on; legend show; xlim([0 2000]);


%% ------------------------- Output Limiter Pass --------------------------
%
% maxItersL = 5;
% peakID = 0;
% sos_limit = [];
% mag_limit_history = [];
% peak_store_limit = [];
% % Compute the residual above the limiter envelope
% limit_resid = H_plot - limit_curve;  % positive values = violating boosts
% mag_limit_current = H_plot;          % start from the applied boost response
%
% figure('Name','Output Limiting','NumberTitle','off');
% semilogx(freq, limit_curve, '--k', 'LineWidth', 1.5, 'DisplayName', 'Limiter Envelope');
% hold on;
% semilogx(freq, H_plot, '-.', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, 'DisplayName', 'Input Filters');
% xlabel('Frequency (Hz)');
% ylabel('Boost (dB)');
% grid on; legend show;
% xlim([0 2000]);
%
% for iter = 1:maxItersL
%
%     if all(limit_resid <= 0)
%         fprintf('No peaks above limiter found. Stopping at iteration %d.\n', iter-1);
%         break;
%     end
%
%
%     lim_peak = findWorstPeak(limit_resid, freq, 0.5, correction_range);
%
%
%     if isempty(lim_peak) || lim_peak.amp < 0.3
%         fprintf('Limiter converged after %d iterations.\n', iter-1);
%         break;
%     end
%
%
%     f0   = lim_peak.f0;
%     gain = -lim_peak.amp;
%     Qval = lim_peak.Q;
%
%     % Design reduction filter
%     sos = designPeakingBiquad(f0, fs, gain, Qval);
%
%     % Apply filter
%     mag_limit_current = applyFilter(mag_limit_current, freq, fs, sos);
%     fprintf('Filter applied at %.1f Hz | Q = %.2f | gain = %.2f dB\n', f0, Qval, gain);
%
%     % Update the residual
%     limit_resid = mag_limit_current - limit_curve;
%
%     % Store filter data
%     peakID = peakID + 1;
%     peak_store_limit(iter,:) = [f0, gain, Qval, peakID];
%     sos_limit = [sos_limit; sos];
%     mag_limit_history(:,iter) = mag_limit_current;
%
%     % --- Plot evolution ---
%     if iter > 1
%         for k = 1:iter-1
%             semilogx(freq, mag_limit_history(:,k), '--', ...
%                 'Color', [0.8 0.8 0.8], 'LineWidth', 0.8, ...
%                 'HandleVisibility','off');
%         end
%     end
%
%     semilogx(freq, mag_limit_current, 'LineWidth', 1.6, ...
%         'Color', [0 0.45 0.9], ...
%         'DisplayName', sprintf('Limited Response %d', iter));
%
%     xlabel('Frequency (Hz)'); ylabel('Boost (dB)');
%     title('Output Limiting Progress');
%     grid on; legend show;
%     drawnow;
% end
%
% peak_store = [peak_store1; peak_store_limit];

maxItersL = 5;
peakID = 0;
sos_limit = [];
mag_limit_history = [];
peak_store_limit = [];

% Compute the residual above the limiter envelope
limit_resid = H_plot - limit_curve;  % positive values = violating boosts
mag_limit_current = H_plot;          % start from the applied boost response

figure('Name','Output Limiting','NumberTitle','off');
semilogx(freq, limit_curve, '--k', 'LineWidth', 1.5, 'DisplayName', 'Limiter Envelope');
hold on;
semilogx(freq, H_plot, '-.', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, 'DisplayName', 'Input Filters');
xlabel('Frequency (Hz)');
ylabel('Boost (dB)');
grid on; legend show;
xlim([0 2000]);

for iter = 1:maxItersL
    if all(limit_resid <= 0)
        fprintf('No peaks above limiter found. Stopping at iteration %d.\n', iter-1);
        break;
    end

    lim_peak = findWorstPeak(limit_resid, freq, 0.5, correction_range);

    if isempty(lim_peak) || lim_peak.amp < 0.3
        fprintf('Limiter converged after %d iterations.\n', iter-1);
        break;
    end

    f0   = lim_peak.f0;
    amp_violation = lim_peak.amp;  % How much we're over the limit

    % === MEASURE THE ACTUAL PEAK WIDTH ===
    % Find where the violation actually exists
    [~, peak_idx] = min(abs(freq - f0));

    % Find the bandwidth where we're actually above the limit
    violation_threshold = 0.1;  % only care about violations > 0.1 dB

    % Search left from peak
    iL = peak_idx;
    while iL > 1 && limit_resid(iL) > violation_threshold
        iL = iL - 1;
    end
    fL_violation = freq(iL);

    % Search right from peak
    iR = peak_idx;
    while iR < length(freq) && limit_resid(iR) > violation_threshold
        iR = iR + 1;
    end
    fR_violation = freq(iR);

    % Calculate the actual violation bandwidth
    BW_violation = fR_violation - fL_violation;
    BW_oct_violation = log2(fR_violation / fL_violation);

    % Calculate Q needed to match this bandwidth
    % For a peaking filter, -3dB bandwidth ≈ f0/Q
    % We want even narrower, so use the violation bandwidth directly
    Q_needed = f0 / BW_violation;
    Q_needed = max(Q_needed, 10);  % minimum Q of 10
    Q_needed = min(Q_needed, 30);  % maximum Q of 30 for stability

    fprintf('\nIter %d Analysis:\n', iter);
    fprintf('  Peak at %.1f Hz, violation %.2f dB\n', f0, amp_violation);
    fprintf('  Violation bandwidth: %.1f-%.1f Hz (%.3f octaves)\n', ...
        fL_violation, fR_violation, BW_oct_violation);
    fprintf('  Required Q: %.2f\n', Q_needed);

    % === DESIGN SURGICAL FILTER ===
    % Use exact gain needed at the peak
    gain_needed = -amp_violation * 1.1;  % 10% extra to ensure we clear the limit

    % Design the filter
    sos = designPeakingBiquad(f0, fs, gain_needed, Q_needed);

    % === VERIFY FILTER RESPONSE ===
    % Calculate what this filter actually does across frequency
    filter_response = applyFilter(zeros(size(freq)), freq, fs, sos);

    % Find -3dB bandwidth of the filter
    peak_filter_response = min(filter_response);  % most negative point
    half_gain = peak_filter_response / 2;

    iL_filter = peak_idx;
    while iL_filter > 1 && filter_response(iL_filter) < half_gain
        iL_filter = iL_filter - 1;
    end

    iR_filter = peak_idx;
    while iR_filter < length(freq) && filter_response(iR_filter) < half_gain
        iR_filter = iR_filter + 1;
    end

    BW_filter = freq(iR_filter) - freq(iL_filter);
    BW_oct_filter = log2(freq(iR_filter) / freq(iL_filter));

    fprintf('  Filter -3dB bandwidth: %.1f-%.1f Hz (%.3f octaves)\n', ...
        freq(iL_filter), freq(iR_filter), BW_oct_filter);

    % === CHECK COLLATERAL DAMAGE BEFORE APPLYING ===
    mag_test = applyFilter(mag_limit_current, freq, fs, sos);

    % Define "safe zone" - slightly wider than violation zone
    safe_margin = 1.5;
    f_safe_lo = fL_violation / safe_margin;
    f_safe_hi = fR_violation * safe_margin;

    % Check impact outside safe zone
    outside_safe = (freq < f_safe_lo) | (freq > f_safe_hi);
    outside_safe = outside_safe & (freq >= correction_range(1)) & (freq <= correction_range(2));

    collateral_change = abs(mag_test(outside_safe) - mag_limit_current(outside_safe));
    max_collateral = max(collateral_change);
    mean_collateral = mean(collateral_change);

    fprintf('  Collateral damage outside %.1f-%.1f Hz:\n', f_safe_lo, f_safe_hi);
    fprintf('    Max: %.3f dB, Mean: %.3f dB\n', max_collateral, mean_collateral);

    % === APPLY FILTER ONLY IF ACCEPTABLE ===
    % If collateral damage is too high relative to the fix, warn user
    if max_collateral > amp_violation * 0.3
        fprintf('  WARNING: High collateral damage (%.3f dB vs %.3f dB violation)\n', ...
            max_collateral, amp_violation);
        fprintf('  Consider if this limiter filter is appropriate.\n');
    end

    % Apply filter
    mag_before_limit = mag_limit_current;
    mag_limit_current = mag_test;

    % Verify we actually reduced the violation
    new_violation_at_peak = mag_limit_current(peak_idx) - limit_curve(peak_idx);
    actual_reduction = (mag_before_limit(peak_idx) - limit_curve(peak_idx)) - new_violation_at_peak;

    fprintf('  Actual peak reduction: %.3f dB (from %.3f to %.3f dB above limit)\n', ...
        actual_reduction, mag_before_limit(peak_idx) - limit_curve(peak_idx), ...
        max(0, new_violation_at_peak));

    % Update the residual
    limit_resid = mag_limit_current - limit_curve;

    % Store filter data
    peakID = peakID + 1;
    peak_store_limit(iter,:) = [f0, gain_needed, Q_needed, peakID];
    sos_limit = [sos_limit; sos];
    mag_limit_history(:,iter) = mag_limit_current;

    % --- DETAILED PLOT ---
    clf;

    % Main plot
    subplot(2,1,1);
    semilogx(freq, limit_curve, '--k', 'LineWidth', 1.5, 'DisplayName', 'Limiter Envelope');
    hold on;

    if iter > 1
        for k = 1:iter-1
            semilogx(freq, mag_limit_history(:,k), '--', ...
                'Color', [0.8 0.8 0.8], 'LineWidth', 0.8, ...
                'HandleVisibility','off');
        end
    end

    semilogx(freq, mag_before_limit, ':', 'Color', [0.6 0.6 0.6], ...
        'LineWidth', 1.2, 'DisplayName', 'Before this filter');
    semilogx(freq, mag_limit_current, 'LineWidth', 1.8, ...
        'Color', [0 0.45 0.9], ...
        'DisplayName', sprintf('After filter %d', iter));

    % Highlight violation zone (red)
    xregion(fL_violation, fR_violation, 'FaceColor', [1 0.8 0.8], 'FaceAlpha', 0.4, ...
        'HandleVisibility', 'off');

    % Highlight safe zone (yellow)
    xregion(f_safe_lo, fL_violation, 'FaceColor', [1 1 0.8], 'FaceAlpha', 0.3, ...
        'HandleVisibility', 'off');
    xregion(fR_violation, f_safe_hi, 'FaceColor', [1 1 0.8], 'FaceAlpha', 0.3, ...
        'HandleVisibility', 'off');

    % Mark the filter location
    xline(f0, 'r--', sprintf('%.0f Hz', f0), 'HandleVisibility', 'off');

    xlabel('Frequency (Hz)'); ylabel('Boost (dB)');
    title(sprintf('Iter %d: f0=%.1f Hz, Q=%.1f, Gain=%.2f dB | Collateral: %.3f dB max', ...
        iter, f0, Q_needed, gain_needed, max_collateral));
    legend('Location', 'best');
    grid on;
    xlim([20 2000]);

    % Filter response plot
    subplot(2,1,2);
    semilogx(freq, filter_response, 'b-', 'LineWidth', 1.5, ...
        'DisplayName', 'This filter response');
    hold on;
    yline(peak_filter_response/2, 'k--', '-3dB', 'DisplayName', '-3dB line');

    % Highlight filter bandwidth
    xregion(freq(iL_filter), freq(iR_filter), 'FaceColor', [0.8 0.9 1], ...
        'FaceAlpha', 0.4, 'DisplayName', 'Filter -3dB BW');

    % Show violation zone for comparison
    xregion(fL_violation, fR_violation, 'FaceColor', [1 0.8 0.8], 'FaceAlpha', 0.3, ...
        'DisplayName', 'Violation zone');

    xline(f0, 'r--', 'HandleVisibility', 'off');

    xlabel('Frequency (Hz)');
    ylabel('Filter Gain (dB)');
    title(sprintf('Filter Response: BW=%.1f Hz (%.3f oct), Q=%.1f', ...
        BW_filter, BW_oct_filter, Q_needed));
    legend('Location', 'best');
    grid on;
    xlim([20 2000]);

    drawnow;
end

fprintf('\n=== Limiter Summary ===\n');
fprintf('Total limiter filters: %d\n', size(sos_limit, 1));
if ~isempty(peak_store_limit)
    fprintf('Average Q: %.2f\n', mean(peak_store_limit(:,3)));
    fprintf('Q range: %.2f - %.2f\n', min(peak_store_limit(:,3)), max(peak_store_limit(:,3)));
    fprintf('Average gain: %.2f dB\n', mean(peak_store_limit(:,2)));
end

peak_store = [peak_store1; peak_store_limit];
%% --------------------- Final Output Plotting ----------------------------

sos_total = [sos_all; sos_limit];
nFilters  = size(sos_total,1);

% expected responce with limiters
mag_expected = mag;
for i = 1:size(sos_total, 1)
    sos_i = sos_total(i, :);
    mag_expected = applyFilter(mag_expected, freq, fs, sos_i);
end


figure('Name','expected responce with limits')
semilogx(freq, L_combined, 'g--','LineWidth',1.5,'DisplayName', 'Target');
hold on;
% semilogx(freq, L_combined+3, '--','DisplayName', 'Target + 3','Color',[0.8 0.5 0.5]);
% semilogx(freq, L_combined-3, '--','DisplayName', 'Target - 3','Color',[0.8 0.5 0.5]);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title(sprintf('Expected Response - %d Filters Applied', size(sos_all, 1)));
legend show;
grid minor;
xlim([20 1000]);
ylim([30 max(mag)+10])

mag_iter = mixedSmooth;

for i = 1:nFilters

    sos_i = sos_total(i,:);
    prev_mag = mag_iter;
    mag_iter = applyFilter(prev_mag, freq, fs, sos_i);

    % Plot current iteration
    if i < nFilters
        semilogx(freq, mag_iter, 'Color',[0.9 0.9 0.9], ...
            'LineWidth',1.0, ...
            'DisplayName',sprintf('After Filter %d', i),'HandleVisibility','off');
    else % for final magnitude responce
        semilogx(freq,PPB,'Color',[1, 0.647, 0],'LineWidth',1,'DisplayName',sprintf('Pre Limiter (Filter %d)',maxIters));
        % Final iteration: bold and blue
        % semilogx(freq, mag_iter,'Color',[0.5, 0, 0.5],'LineWidth',2.0,'DisplayName',sprintf('Final Output (Filter %d)', i));
    end

    drawnow;
end

semilogx(freq, mag_expected,'r', 'LineWidth',1.8,'DisplayName', 'Mag with applied filters');
semilogx(freq, mag,'k','LineWidth', 1, 'DisplayName', 'REW Input');


%% ------------------------ Export Filters --------------------------------

% Prepare to export filter values
f0s = peak_store(:, 1);
Qvals = peak_store(:, 3);
gains = peak_store(:, 2);

Nf = length(f0s);
if length(Qvals) ~= Nf || length(gains) ~= Nf
    error('f0s, Qvals and gains must be the same length.');
end

outFilename = 'room_eq.txt';

% Open file
fid = fopen(outFilename, 'w');
if fid == -1
    error('Could not open %s for writing.', outFilename);
end

fmt = '%.2f %.2f %.2f\n';


for k = 1:Nf
    fprintf(fid, fmt, round(f0s(k),2,"decimals"), ...
        round(Qvals(k),2,"decimals"), ...
        round(gains(k),2,"decimals") ...
        );
end

% pad last filters:

for k = (Nf+1):20

    fprintf(fid,fmt,1000,1,0);

end

fclose(fid);
fprintf('Wrote %d filters to: %s\n', Nf, outFilename);



%% ---------------------------- Functions ---------------------------------
function sos = designPeakingBiquad(f0, fs, dBgain, Q)
A = 10^(dBgain/40);
w0 = 2*pi*f0/fs;
alpha = sin(w0)/(2*Q);

b0 = 1 + alpha*A;
b1 = -2*cos(w0);
b2 = 1 - alpha*A;

a0 = 1 + alpha/A;
a1 = -2*cos(w0);
a2 = 1 - alpha/A;

sos = [b0 b1 b2 a0 a1 a2];
end

function [pk, lc] = s_findpeaks(x, thresh, minDist)
N = length(x);
pk = [];
lc = [];

% find peaks peaks
for i = 2:N-1  % skip first and last sample
    if x(i) > thresh && x(i) > x(i-1) && x(i) > x(i+1)
        pk(end+1) = x(i); % store peak value
        lc(end+1) = i; % store peak index
    end
end

% add minimum distance between peaks
if ~isempty(lc)
    keep = true(size(lc));
    last = lc(1);  % last accepted peak index

    for k = 2:length(lc) % for all peaks
        if lc(k) - last < minDist  % find if difference in freq is less than minDist
            if pk(k) > pk(k-1) % find out which peak is bigger
                keep(k-1) = false; % remove the smaller previous peak
                last = lc(k);
            else
                keep(k) = false;   % remove current smaller peak
            end
        else
            last = lc(k); % update last accepted index
        end
    end

    pk = pk(keep); % keep only selected peaks
    lc = lc(keep); % keep only selected peak frequency
end
end


%
% function [freq, mag] = processWavToRTA(wavFilePath, varargin)
% % processWavToRTA  Analyze room frequency response from recorded test signal
% %
% % This processes a microphone recording of pink noise or sweep played in a room,
% % similar to REW's RTA (Real-Time Analyzer) function.
% %
% % INPUTS:
% %   wavFilePath - path to WAV file of recorded test signal
% %
% % OPTIONAL NAME-VALUE PAIRS:
% %   'ReferenceFile' - path to reference signal WAV (for transfer function)
% %   'WindowSize'    - FFT window size in seconds (default: 0.5)
% %   'Overlap'       - overlap fraction (default: 0.75)
% %   'Smoothing'     - fractional octave smoothing (default: 1/48)
% %   'FreqRange'     - [fmin fmax] frequency range (default: [20 2000])
% %
% % OUTPUTS:
% %   freq - frequency vector (Hz)
% %   mag  - magnitude vector (dB)
% %
% % EXAMPLES:
% %   % Basic RTA from pink noise recording
% %   [freq, mag] = processWavToRTA('room_recording.wav');
% %
% %   % With reference signal (for transfer function)
% %   [freq, mag] = processWavToRTA('room_recording.wav', 'ReferenceFile', 'pink_noise_source.wav');
%
% % Parse inputs
% p = inputParser;
% addParameter(p, 'ReferenceFile', '');
% addParameter(p, 'WindowSize', 0.5);      % 0.5 seconds
% addParameter(p, 'Overlap', 0.75);        % 75% overlap
% addParameter(p, 'Smoothing', 1/48);      % 1/48 octave
% addParameter(p, 'FreqRange', [20 2000]);
% parse(p, varargin{:});
%
% opts = p.Results;
%
% % --- Load recorded signal ---
% [x, fs] = audioread(wavFilePath);
% x = x(:,1);  % mono
% x = x(:);    % column vector
%
% fprintf('Processing: %s\n', wavFilePath);
% fprintf('Sample rate: %d Hz, Duration: %.2f s\n', fs, length(x)/fs);
%
% % --- Parameters based on window size ---
% windowTime = opts.WindowSize;
% Nfft = 2^nextpow2(windowTime * fs);  % Next power of 2
% overlap_frac = opts.Overlap;
% hop = max(1, round(Nfft * (1 - overlap_frac)));
%
% fprintf('FFT size: %d points (%.3f s window)\n', Nfft, Nfft/fs);
% fprintf('Overlap: %.1f%% (%d samples hop)\n', overlap_frac*100, hop);
%
% % --- Hann window ---
% n = (0:Nfft-1)';
% win = 0.5 - 0.5*cos(2*pi*n/(Nfft-1));
%
% % --- Frame the signal ---
% L = length(x);
% nFrames = floor((L - Nfft) / hop) + 1;
%
% if nFrames < 1
%     error('Signal too short for the specified window size');
% end
%
% fprintf('Processing %d frames...\n', nFrames);
%
% % --- Compute average spectrum ---
% half = Nfft/2 + 1;
% Psum = zeros(half, 1);
%
% for k = 1:nFrames
%     startIdx = (k-1)*hop + 1;
%     endIdx = startIdx + Nfft - 1;
%
%     if endIdx > L
%         break;
%     end
%
%     frame = x(startIdx:endIdx) .* win;
%     X = fft(frame);
%
%     % Power spectrum
%     mag2 = abs(X).^2;
%     P = mag2(1:half);
%
%     % Double non-DC/non-Nyquist bins (single-sided spectrum)
%     if half > 2
%         P(2:end-1) = 2 * P(2:end-1);
%     end
%
%     Psum = Psum + P;
% end
%
% % --- Average and convert to dB ---
% Pavg = Psum / nFrames;
% MAG = sqrt(Pavg);
% MAGdB = 20*log10(MAG + 1e-12);
%
% % --- Frequency vector ---
% freq = (0:half-1)' * (fs/Nfft);
%
% % Remove DC (0 Hz)
% freq = freq(2:end);
% MAGdB = MAGdB(2:end);
%
% % --- Apply frequency range ---
% fRange = opts.FreqRange;
% valid_idx = freq >= fRange(1) & freq <= fRange(2);
% freq = freq(valid_idx);
% MAGdB = MAGdB(valid_idx);
%
% % --- Fractional octave smoothing ---
% N_oct = 1 / opts.Smoothing;  % e.g., 1/48 → 48
% sigma = 1 / (2 * N_oct);
%
% lf = log2(freq);
% MAGdB_smooth = zeros(size(MAGdB));
%
% fprintf('Applying 1/%.0f octave smoothing...\n', N_oct);
%
% for i = 1:length(freq)
%     d = lf - lf(i);
%     w = exp(-(d.^2) / (2*sigma^2));
%     w = w / sum(w);
%     MAGdB_smooth(i) = sum(w .* MAGdB);
% end
%
% mag = MAGdB_smooth;
%
% % --- Clean up NaN/Inf ---
% valid = isfinite(freq) & isfinite(mag);
% if ~all(valid)
%     warning('Removed %d non-finite values', sum(~valid));
% end
%
% freq = freq(valid);
% mag = mag(valid);
%
% % --- Optional: Compute transfer function with reference ---
% if ~isempty(opts.ReferenceFile)
%     fprintf('Computing transfer function with reference: %s\n', opts.ReferenceFile);
%     [x_ref, fs_ref] = audioread(opts.ReferenceFile);
%
%     if fs_ref ~= fs
%         error('Reference file sample rate (%d) does not match recording (%d)', fs_ref, fs);
%     end
%
%     x_ref = x_ref(:,1);
%     x_ref = x_ref(:);
%
%     % Process reference the same way
%     Psum_ref = zeros(half, 1);
%     nFrames_ref = floor((length(x_ref) - Nfft) / hop) + 1;
%
%     for k = 1:nFrames_ref
%         startIdx = (k-1)*hop + 1;
%         endIdx = startIdx + Nfft - 1;
%
%         if endIdx > length(x_ref)
%             break;
%         end
%
%         frame = x_ref(startIdx:endIdx) .* win;
%         X = fft(frame);
%         mag2 = abs(X).^2;
%         P = mag2(1:half);
%
%         if half > 2
%             P(2:end-1) = 2 * P(2:end-1);
%         end
%
%         Psum_ref = Psum_ref + P;
%     end
%
%     Pavg_ref = Psum_ref / nFrames_ref;
%     MAG_ref = sqrt(Pavg_ref);
%     MAGdB_ref = 20*log10(MAG_ref + 1e-12);
%
%     % Remove DC
%     MAGdB_ref = MAGdB_ref(2:end);
%     MAGdB_ref = MAGdB_ref(valid_idx);
%
%     % Smooth reference
%     MAGdB_ref_smooth = zeros(size(MAGdB_ref));
%     for i = 1:length(freq)
%         d = lf - lf(i);
%         w = exp(-(d.^2) / (2*sigma^2));
%         w = w / sum(w);
%         MAGdB_ref_smooth(i) = sum(w .* MAGdB_ref);
%     end
%
%     % Transfer function = measured / reference
%     mag = mag - MAGdB_ref_smooth;
%     fprintf('Transfer function computed (room response relative to source)\n');
% end
%
% fprintf('RTA complete: %d frequency points from %.1f to %.1f Hz\n', ...
%     length(freq), freq(1), freq(end));
%
% end

function [freq, mag] = processWavToRTA(wavFilePath)
% processWavToRTA  REW-style RTA analysis with averaging over entire recording
%
% Processes a WAV file recording using REW RTA settings:
% - 32k FFT (32768 points)
% - Hann window
% - 93.75% overlap (6.25% hop = 2048 samples at 44.1kHz)
% - 1/48 octave smoothing
% - Averages all frames across the entire recording
%
% INPUTS:
%   wavFilePath - path to WAV file recording
%
% OUTPUTS:
%   freq - frequency vector (Hz)
%   mag  - magnitude vector (dB)
%
% EXAMPLE:
%   [freq, mag] = processWavToRTA('My recording 10.wav');

% --- Load WAV file ---
[x, fs] = audioread(wavFilePath);
x = x(:,1);  % mono
x = x(:);    % column vector

fprintf('Processing: %s\n', wavFilePath);
fprintf('Sample rate: %d Hz\n', fs);
fprintf('Duration: %.2f seconds\n', length(x)/fs);

% --- REW RTA Parameters ---
Nfft = 32768;                    % 32k FFT
overlap_frac = 0.9375;           % 93.75% overlap
hop = round(Nfft * (1 - overlap_frac));  % 2048 samples

fprintf('FFT size: %d points (%.3f second window)\n', Nfft, Nfft/fs);
fprintf('Overlap: %.1f%% (%d samples hop)\n', overlap_frac*100, hop);

% --- Hann window ---
n = (0:Nfft-1)';
win = 0.5 - 0.5*cos(2*pi*n/(Nfft-1));

% --- Calculate number of frames ---
L = length(x);
nFrames = floor((L - Nfft) / hop) + 1;

if nFrames < 1
    error('Recording too short for 32k FFT window');
end

fprintf('Total frames to average: %d\n', nFrames);

% --- Process all frames and accumulate ---
half = Nfft/2 + 1;
Psum = zeros(half, 1);

for k = 1:nFrames
    startIdx = (k-1)*hop + 1;
    endIdx = startIdx + Nfft - 1;

    if endIdx > L
        break;
    end

    % Extract frame and apply window
    frame = x(startIdx:endIdx) .* win;

    % FFT
    X = fft(frame);

    % Power spectrum (magnitude squared)
    mag2 = abs(X).^2;
    P = mag2(1:half);

    % Double non-DC/non-Nyquist bins for single-sided spectrum
    if half > 2
        P(2:end-1) = 2 * P(2:end-1);
    end

    % Accumulate
    Psum = Psum + P;
end

% --- Average power across all frames ---
Pavg = Psum / nFrames;

% --- Convert to magnitude in dB ---
MAG = sqrt(Pavg);                    % RMS magnitude
MAGdB = 20*log10(MAG + 1e-12);       % Convert to dB

% --- Create frequency vector ---
freq = (0:half-1)' * (fs/Nfft);

% --- Remove DC component (0 Hz) ---
freq = freq(2:end);
MAGdB = MAGdB(2:end);

% --- Restrict to 20-2000 Hz range ---
valid_idx = freq >= 20 & freq <= 2000;
freq = freq(valid_idx);
MAGdB = MAGdB(valid_idx);

% --- Apply 1/48 octave smoothing (REW style) ---
N_oct = 48;                          % 1/48 octave
sigma = 1/(2*N_oct);                 % Gaussian width in octaves

lf = log2(freq);                     % Log frequency scale
MAGdB_smooth = zeros(size(MAGdB));

fprintf('Applying 1/48 octave smoothing...\n');

for i = 1:length(freq)
    d = lf - lf(i);                  % Distance in octaves
    w = exp(-(d.^2) / (2*sigma^2));  % Gaussian weights
    w = w / sum(w);                  % Normalize
    MAGdB_smooth(i) = sum(w .* MAGdB);
end

mag = MAGdB_smooth;

% --- Remove any NaN or Inf values ---
valid = isfinite(freq) & isfinite(mag);

if ~all(valid)
    warning('Removed %d non-finite values from output', sum(~valid));
end

freq = freq(valid);
mag = mag(valid);

% --- Summary ---
fprintf('\n=== Processing Complete ===\n');
fprintf('Averaged %d frames over %.2f seconds\n', nFrames, length(x)/fs);
fprintf('Output: %d frequency points from %.1f to %.1f Hz\n', ...
    length(freq), freq(1), freq(end));
fprintf('Frequency resolution: %.2f Hz\n', fs/Nfft);

end





% function [freq, mag] = processWavToFreqMag(wavFilePath)
% % processWavToFreqMag  Convert WAV file to frequency-magnitude response
% %
% % INPUTS:
% %   wavFilePath - path to WAV file (string or char)
% %
% % OUTPUTS:
% %   freq - frequency vector (Hz), cleaned of NaN/Inf
% %   mag  - magnitude vector (dB), cleaned of NaN/Inf
% %
% % EXAMPLE:
% %   [freq, mag] = processWavToFreqMag('white_noise_263s_Matlab10.wav');
%
% % --- WAV file to array ---
% [x, fs] = audioread(wavFilePath);
% x = x(:,1);            % mono
% x = x(:);              % column vector
%
% % --- Parameters ---
% Nfft = 32768;          % 32k FFT
% overlap_frac = 0.93;   % 93% overlap
% hop = max(1, round(Nfft * (1-overlap_frac)));   % hop = 7%
%
% % --- Manual Hann window ---
% n = (0:Nfft-1)';
% win = 0.5 - 0.5*cos(2*pi*n/(Nfft-1));   % manual Hann
%
% % --- Manual framing ---
% L = length(x);
% nFrames = ceil((L - Nfft) / hop) + 1;
% padNeeded = (nFrames-1)*hop + Nfft - L;
% xpad = [x; zeros(padNeeded,1)];
%
% % --- Manual FFT loop ---
% half = Nfft/2 + 1;
% Psum = zeros(half,1);
%
% for k = 1:nFrames
%     idx = (k-1)*hop + (1:Nfft);
%     frame = xpad(idx) .* win;
%     X = fft(frame);               % full FFT
%     mag2 = abs(X).^2;             % magnitude^2 = power
%     P = mag2(1:half);             % single-sided power
%
%     % Double non-DC/non-Nyquist to restore lost negative frequencies
%     if half > 2
%         P(2:end-1) = 2 * P(2:end-1);
%     end
%
%     Psum = Psum + P;
% end
%
% % --- Average power ---
% Pavg = Psum / nFrames;
% MAG = sqrt(Pavg);                     % linear magnitude
% MAGdB = 20*log10(MAG + 1e-12);        % in dB
% freq = (0:half-1)' * (fs/Nfft);       % frequency axis
%
% % --- Frequency limit to speed up calc ---
% freqlim = 2000;
% valid_idx = freq <= freqlim;
% freq = freq(valid_idx);
% MAGdB = MAGdB(valid_idx);
%
% % --- 1/48-Oct smoothing (manual, Gaussian in log-frequency) ---
% N_oct = 48;                           % 1/48 octave
% sigma = 1/(2*N_oct);                  % Gaussian width in octaves
% lf = log2(freq);
% lf(1) = lf(2) - (lf(3)-lf(2));        % fix DC (avoid -inf)
%
% MAGdB_smooth = zeros(size(MAGdB));
% for i = 1:length(freq)
%     d = lf - lf(i);                   % difference in octaves
%     w = exp(-(d.^2) / (2*sigma^2));   % Gaussian weights
%     w = w / sum(w);                   % normalize
%     MAGdB_smooth(i) = sum(w .* MAGdB);
% end
%
% % --- Use smoothed magnitude ---
% mag = MAGdB_smooth;
%
% % --- Remove NaN and Inf values ---
% valid = isfinite(freq) & isfinite(mag);
%
% if ~all(valid)
%     warning('Removed %d non-finite values from output', sum(~valid));
% end
%
% freq = freq(valid);
% mag = mag(valid);
%
% % --- Final check ---
% if any(~isfinite(freq)) || any(~isfinite(mag))
%     error('Output still contains NaN or Inf values after cleaning');
% end
%
% fprintf('Successfully processed: %s\n', wavFilePath);
% fprintf('Output: %d frequency points from %.1f to %.1f Hz\n', ...
%     length(freq), freq(1), freq(end));
%
% end

function peak = findWorstPeak(resid, freq, thresh, freqRange)

% Default range if not specified
if nargin < 4
    freqRange = [20 500];
end

if nargin < 5
    maxQ = 15;  % maximum Q
end

% Restrict search range
idx_range = (freq >= freqRange(1) & freq <= freqRange(2));
freq_sub  = freq(idx_range);
resid_sub = resid(idx_range);

% Find peaks in restricted range

% [pk, lc] = findpeaks(resid_sub, "MinPeakHeight", thresh, "MinPeakDistance", 10);
[pk, lc] = s_findpeaks(resid_sub, thresh, 10);
if isempty(pk)
    peak = [];
    return;
end

% find largest peak
[~, idxMax] = max(pk);
amp = pk(idxMax);
idx = lc(idxMax);
f0  = freq_sub(idx);

% find Q by half-height bandwidth
halfHeight = amp/2;
iL = idx; iR = idx;
while iL > 1 && resid_sub(iL) > halfHeight, iL = iL - 1; end
while iR < length(freq_sub) && resid_sub(iR) > halfHeight, iR = iR + 1; end

% if iL < idx && iR > idx
BW = freq_sub(iR) - freq_sub(iL);
Q  = f0 / BW;
Q = min(Q, maxQ);
% else
%     Q = 6; % fallback
% end

peak = struct("f0", f0, "amp", amp, "Q", Q);


end

function mag_next = applyFilter(mag_current, freq, fs, sos)

% Extract normalized biquad coefficients
b = sos(1:3)/sos(4);
a = sos(4:6)/sos(4);

% Convert frequency in Hz to radians
w = 2*pi*freq./fs;

% Compute H(e^jw) manually (complex response)
ejw = exp(-1j*w);
num = b(1) + b(2)*ejw + b(3)*ejw.^2;
den = a(1) + a(2)*ejw + a(3)*ejw.^2;
H = num ./ den;

% Apply magnitude response
mag_lin = 10.^(mag_current/20) .* abs(H);
mag_next = 20*log10(mag_lin);

end

function plotIterationProgress(iter, freq, mag_orig, mag_current, target, sos_all)
% Faster version of the original plotIterationProgress()
% Keeps identical visuals but updates in-place with semilogx scale.

persistent hFig hAx hLines

if isempty(hFig) || ~isvalid(hFig)
    % === Initialize once (first iteration) ===
    hFig = figure('Name','Iteration Progress','NumberTitle','off');
    clf(hFig);

    % --- Top subplot (Response vs Target) ---
    hAx(1) = subplot(2,1,1,'Parent',hFig);
    hold(hAx(1), 'on');
    hLines.orig = semilogx(hAx(1), freq, mag_orig, ...
        'DisplayName','Original', 'LineWidth',1);
    hLines.curr = semilogx(hAx(1), freq, mag_current, ...
        'DisplayName','Current', 'LineWidth',1.5);
    hLines.targ = semilogx(hAx(1), freq, target, ...
        'g--','DisplayName','Target','LineWidth',1.2);
    grid(hAx(1), 'on');
    xlim(hAx(1), [20 1000]);
    xlabel(hAx(1), 'Frequency (Hz)');
    ylabel(hAx(1), 'Magnitude (dB)');
    title(hAx(1), sprintf('Iteration %d - Response vs Target', iter));
    legend(hAx(1), 'show');

    % Force semilog scale (even if reused)
    set(hAx(1), 'XScale', 'log');

    % --- Bottom subplot (Applied Reduction) ---
    hAx(2) = subplot(2,1,2,'Parent',hFig);
    hold(hAx(2), 'on');
    hLines.reduction = semilogx(hAx(2), freq, (mag_orig - mag_current)*-1, ...
        'm','DisplayName','Reduction','LineWidth',1.2);
    grid(hAx(2), 'on');
    xlim(hAx(2), [20 1000]);
    xlabel(hAx(2), 'Frequency (Hz)');
    ylabel(hAx(2), 'Reduced level (dB)');
    title(hAx(2), 'Applied Reduction at Peaks');
    legend(hAx(2), 'show');
    set(hAx(2), 'XScale', 'log'); % enforce semilogx

else
    % === Fast update path ===
    set(hLines.curr, 'YData', mag_current);
    set(hLines.targ, 'YData', target);
    set(hLines.reduction, 'YData', (mag_orig - mag_current)*-1);
    title(hAx(1), sprintf('Iteration %d - Response vs Target', iter));

    % Ensure semilog scaling stays correct
    set(hAx(1), 'XScale', 'log');
    set(hAx(2), 'XScale', 'log');
end

drawnow limitrate nocallbacks;
pause(0.3)
end


function limitConfig = generateLimitBands(f_start, f_end, n_bands, varargin)
% generateLimitBands  Automatically create frequency bands with smooth limit transitions
%
% INPUTS:
%   f_start  - start frequency (Hz), e.g., 20
%   f_end    - end frequency (Hz), e.g., 500
%   n_bands  - number of bands to create
%
% OPTIONAL NAME-VALUE PAIRS:
%   'spacing'       - 'log' (default) or 'linear'
%   'boost_start'   - max boost at f_start (dB), default 2.5
%   'boost_end'     - max boost at f_end (dB), default 7.0
%   'cut_start'     - max cut at f_start (dB), default 4.0
%   'cut_end'       - max cut at f_end (dB), default 8.5
%   'Q_min_start'   - min Q at f_start, default 0.5
%   'Q_min_end'     - min Q at f_end, default 2.5
%   'Q_max_start'   - max Q at f_start, default 2.0
%   'Q_max_end'     - max Q at f_end, default 15.0
%   'energy_start'  - energy limit at f_start, default 8
%   'energy_end'    - energy limit at f_end, default 70
%   'transition'    - 'linear', 'exponential', or 'sigmoid', default 'exponential'
%   'f_high_end'    - upper frequency limit for final band (default 2000)
%
% OUTPUT:
%   limitConfig - struct with all band definitions
%
% EXAMPLES:
%   % Conservative 10 bands from 20-500 Hz:
%   cfg = generateLimitBands(20, 500, 10);
%
%   % Aggressive 15 bands with custom limits:
%   cfg = generateLimitBands(20, 500, 15, 'boost_start', 3.5, 'boost_end', 8.0);
%
%   % Linear spacing instead of logarithmic:
%   cfg = generateLimitBands(20, 500, 8, 'spacing', 'linear');

% Parse optional inputs
p = inputParser;
addParameter(p, 'spacing', 'log');
addParameter(p, 'boost_start', 2.5);
addParameter(p, 'boost_end', 7.0);
addParameter(p, 'cut_start', 4.0);
addParameter(p, 'cut_end', 8.5);
addParameter(p, 'Q_min_start', 0.1);
addParameter(p, 'Q_min_end', 0.1);
addParameter(p, 'Q_max_start', 15);
addParameter(p, 'Q_max_end', 15);
addParameter(p, 'energy_start', 8);
addParameter(p, 'energy_end', 70);
addParameter(p, 'transition', 'exponential');
addParameter(p, 'f_high_end', 1000);
parse(p, varargin{:});

opts = p.Results;

% Generate band edges
if strcmp(opts.spacing, 'log')
    % Logarithmic spacing (better for audio)
    band_edges = logspace(log10(f_start), log10(f_end), n_bands + 1);
else
    % Linear spacing
    band_edges = linspace(f_start, f_end, n_bands + 1);
end

% Create band array [f_min, f_max]
limitConfig.bands = zeros(n_bands + 1, 2);
for i = 1:n_bands
    limitConfig.bands(i, :) = [band_edges(i), band_edges(i+1)];
end
% Add high-frequency band
limitConfig.bands(n_bands + 1, :) = [f_end, opts.f_high_end];

% Generate normalized positions for interpolation (0 to 1)
norm_pos = linspace(0, 1, n_bands + 1);

% Apply transition function
switch opts.transition
    case 'linear'
        t = norm_pos;
    case 'exponential'
        % Exponential curve (more gradual at start, steeper at end)
        t = norm_pos.^1.5;
    case 'sigmoid'
        % S-curve (gradual at both ends, steep in middle)
        t = 1 ./ (1 + exp(-10*(norm_pos - 0.5)));
        t = (t - min(t)) / (max(t) - min(t));  % normalize to 0-1
    otherwise
        t = norm_pos;
end

% Interpolate all parameters
limitConfig.max_boost = opts.boost_start + t' * (opts.boost_end - opts.boost_start);
limitConfig.max_cut = opts.cut_start + t' * (opts.cut_end - opts.cut_start);
limitConfig.Q_min = opts.Q_min_start + t' * (opts.Q_min_end - opts.Q_min_start);
limitConfig.Q_max = opts.Q_max_start + t' * (opts.Q_max_end - opts.Q_max_start);
limitConfig.energy_limit = opts.energy_start + t' * (opts.energy_end - opts.energy_start);

% Store configuration details
limitConfig.description = sprintf('%d bands, %.0f-%.0f Hz, %s spacing, %s transition', ...
    n_bands, f_start, f_end, opts.spacing, opts.transition);

end


function [limited_gain, limited_Q, band_idx] = applyFrequencyLimits(f0, raw_gain, raw_Q, limitConfig)
% applyFrequencyLimits  Apply frequency-dependent limits from configuration
%
% INPUTS:
%   f0          - center frequency (Hz)
%   raw_gain    - desired gain (dB)
%   raw_Q       - desired Q factor
%   limitConfig - configuration from generateLimitBands()
%
% OUTPUTS:
%   limited_gain - gain after limiting (dB)
%   limited_Q    - Q after limiting
%   band_idx     - which band was used

% Find which band this frequency belongs to
band_idx = 0;
for i = 1:size(limitConfig.bands, 1)
    if f0 >= limitConfig.bands(i, 1) && f0 <= limitConfig.bands(i, 2)
        band_idx = i;
        break;
    end
end

% If frequency is outside all bands, use last band
if band_idx == 0
    band_idx = size(limitConfig.bands, 1);
end

% Get limits for this band
max_boost = limitConfig.max_boost(band_idx);
max_cut = limitConfig.max_cut(band_idx);
Q_min = limitConfig.Q_min(band_idx);
Q_max = limitConfig.Q_max(band_idx);
energy_limit = limitConfig.energy_limit(band_idx);

% Apply gain limits
if raw_gain > 0
    limited_gain = min(raw_gain, max_boost);
else
    limited_gain = max(raw_gain, -max_cut);
end

% Apply Q limits
limited_Q = max(Q_min, min(raw_Q, Q_max));

% Apply energy limit
gain_Q_product = abs(limited_gain) * limited_Q;
if gain_Q_product > energy_limit
    limited_Q = energy_limit / abs(limited_gain);
    limited_Q = max(Q_min, min(limited_Q, Q_max));
end
end


% function resid = computeResidual(mag_current, target)
% resid = mag_current - target(:);
% end

% function dip = findWorstDip(resid, freq, thresh, fRange)
% % findWorstDip  Detects the most significant dip region (below target)
% %               and computes its center frequency, amplitude, width, and Q.
% %
% % INPUTS:
% %   resid   - residual (measured - target) in dB
% %   freq    - frequency vector (Hz)
% %   thresh  - minimum dip depth to consider (in dB, negative value)
% %   fRange  - [f_min f_max] frequency range to search (Hz)
% %
% % OUTPUT:
% %   dip     - struct with fields:
% %               f0          - center frequency (Hz)
% %               amp         - dip amplitude (dB, negative)
% %               width       - bandwidth (Hz)
% %               Q           - adaptive quality factor
% %               importance  - combined importance metric
%
% % --- Limit search to given frequency range ---
% idx_range = (freq >= fRange(1) & freq <= fRange(2));
% resid_sub = resid(idx_range);
% freq_sub  = freq(idx_range);
%
% % --- Identify regions below threshold (potential dips) ---
% isDip = resid_sub < thresh;
% d = diff([0; isDip; 0]);
% startIdx = find(d == 1);
% endIdx   = find(d == -1) - 1;
%
% if isempty(startIdx)
%     dip = [];
%     return;
% end
%
% % --- Analyze each dip region ---
% nDips = numel(startIdx);
% dipMeans   = zeros(nDips,1);
% dipCenters = zeros(nDips,1);
% dipWidths  = zeros(nDips,1);
% dipImportance = zeros(nDips,1);
%
% for k = 1:nDips
%     idx = startIdx(k):endIdx(k);
%     dipMeans(k)   = mean(resid_sub(idx));              % average depth (dB)
%     dipCenters(k) = mean(freq_sub(idx));               % center freq (Hz)
%     dipWidths(k)  = freq_sub(endIdx(k)) - freq_sub(startIdx(k)); % bandwidth (Hz)
%
%     dipDepth = abs(min(resid_sub(idx)));
%
%     f1 = freq_sub(startIdx(k));
%     f2 = freq_sub(endIdx(k));
%     width_oct = log2(f2 / f1);
%
%     dipImportance(k) =   width_oct/dipDepth;
%
%     %
%     % dipImportance(k) = abs(dipMeans(k)) / max(dipWidths(k), 1);  % importance metric
%
%
%
%
% end
%
% % --- Pick most significant dip ---
% [~, idxBest] = max(dipImportance);
% f0     = dipCenters(idxBest);
% amp    = dipMeans(idxBest);
% width  = dipWidths(idxBest);
%
% % % --- Adaptive Q calculation ---
% Q_base = f0 / max(width, 1e-6);  % avoid division by zero
%
% % % Frequency weighting:
% % if f0 < 150
% %     Q = Q_base * (f0 / 150)^4;     % wider for low frequencies
% % elseif f0 > 2000
% %     Q = Q_base * (f0 / 2000)^0.3;    % slightly narrower at high end
% % else
% Q = Q_base;                      % midband unchanged
% % end
%
% % % Clamp to safe range
% % Q = max(1, min(10, Q));
%
% % --- Package results ---
% dip = struct( ...
%     "f0", f0, ...
%     "amp", amp, ...
%     "width", width, ...
%     "Q", Q, ...
%     "importance", dipImportance(idxBest) ...
%     );
% end
