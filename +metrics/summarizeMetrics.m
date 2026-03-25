function summary = summarizeMetrics(frame_metrics, frame_time)
% SUMMARIZEMETRICS  Aggregate per-frame metrics into cross-frame statistics.
%
%   SUMMARY = SUMMARIZEMETRICS(FRAME_METRICS, FRAME_TIME) takes the
%   struct array of per-frame metrics produced by computeFrameMetrics
%   and computes summary statistics (mean, std, median, min, max)
%   across all processed frames. Additionally computes derived
%   cross-frame metrics that characterize temporal stability.
%
%   Inputs:
%     frame_metrics - [1 x F] or [F x 1] struct array from
%                     computeFrameMetrics, one entry per frame.
%                     Empty entries (skipped frames) are excluded.
%     frame_time    - [F x 1] Processing time per frame [s].
%                     Used for computational cost reporting.
%
%   Outputs:
%     summary - Struct with fields:
%
%       --- Per-metric aggregates ---
%       For each field in frame_metrics (e.g., 'srr'), the summary
%       contains:
%         .srr_mean, .srr_std, .srr_median, .srr_min, .srr_max
%       This pattern repeats for all numeric fields.
%
%       --- Derived cross-frame metrics ---
%       .contrast_stability  - Coefficient of variation (std/mean)
%                              of contrast_iqr across frames. Low
%                              values indicate the accumulator
%                              maintains consistent dynamic range
%                              regardless of event rate. Key metric
%                              for evaluating divisive normalization.
%       .temporal_ssim_stability - Std of temporal_ssim across
%                              frames. Low values indicate uniformly
%                              smooth temporal transitions.
%       .srr_stability       - Std of srr across frames. For an
%                              auto-thresholding filter (coherence),
%                              this should be lower than for
%                              fixed-threshold filters.
%       .n_valid_frames      - Number of frames with valid metrics
%                              (non-empty, non-NaN srr).
%
%       --- Computational cost ---
%       .time_mean           - Mean processing time per frame [s].
%       .time_std            - Std of processing time per frame.
%       .time_total          - Total processing time [s].
%
%       --- Raw vectors (for plotting) ---
%       .srr_vec             - [F x 1] SRR per frame.
%       .temporal_ssim_vec   - [F x 1] SSIM per frame.
%       .contrast_iqr_vec    - [F x 1] Contrast IQR per frame.
%       .edge_sharpness_vec  - [F x 1] Edge sharpness per frame.
%       .noise_floor_vec     - [F x 1] Noise floor per frame.
%       .n_passed_vec        - [F x 1] Events passed per frame.
%       .n_total_vec         - [F x 1] Total events per frame.
%       .polarity_balance_vec - [F x 1] Polarity balance per frame.
%       .n_clusters_vec      - [F x 1] Cluster count per frame.
%       .mean_gradient_vec   - [F x 1] Mean gradient per frame.
%       .surface_fill_vec    - [F x 1] Surface fill per frame.
%
%   Notes:
%     - NaN values are excluded from all aggregate computations
%       using nanmean, nanstd, etc. This handles the first frame
%       (where temporal_ssim is NaN) and any skipped frames.
%     - The raw vectors are stored for downstream plotting (time
%       series, histograms, condition comparisons).
%     - This function does not produce any plots. Use the raw
%       vectors with your own plotting code or a dedicated
%       plotting function.
%
%   Example:
%     % After the main processing loop:
%     summary = metrics.summarizeMetrics(frame_metrics, frame_time);
%     fprintf('Mean SRR: %.3f +/- %.3f\n', ...
%         summary.srr_mean, summary.srr_std);
%     fprintf('Contrast stability (CV): %.4f\n', ...
%         summary.contrast_stability);
%
%   See also: metrics.computeFrameMetrics

    % ================================================================
    % 0. Filter out empty/invalid frames
    % ================================================================
    valid_mask = ~arrayfun(@(s) isempty(fieldnames(s)) || ...
        (isfield(s, 'n_total') && s.n_total == 0), frame_metrics);

    fm = frame_metrics(valid_mask);
    n_valid = numel(fm);

    summary.n_valid_frames = n_valid;

    if n_valid == 0
        warning('summarizeMetrics:noValidFrames', ...
            'No valid frames found. Returning empty summary.');
        return;
    end

    % ================================================================
    % 1. Extract per-frame vectors
    % ================================================================
    fields_to_extract = { ...
        'srr', 'pixel_retention', 'polarity_balance', ...
        'n_passed', 'n_total', 'n_active_pixels', 'n_retained_pixels', ...
        'edge_sharpness', 'mean_gradient', 'contrast_iqr', ...
        'noise_floor', 'surface_fill', ...
        'temporal_ssim', 'temporal_mae', ...
        'n_clusters', 'mean_cluster_size', 'largest_cluster_frac'};

    for i = 1:numel(fields_to_extract)
        fname = fields_to_extract{i};
        if isfield(fm, fname)
            vec = [fm.(fname)]';
        else
            vec = NaN(n_valid, 1);
        end
        summary.([fname '_vec']) = vec;
    end

    % ================================================================
    % 2. Compute aggregate statistics for each field
    % ================================================================
    agg_fields = { ...
        'srr', 'pixel_retention', 'polarity_balance', ...
        'edge_sharpness', 'mean_gradient', 'contrast_iqr', ...
        'noise_floor', 'surface_fill', ...
        'temporal_ssim', 'temporal_mae', ...
        'n_clusters', 'mean_cluster_size', 'largest_cluster_frac'};

    for i = 1:numel(agg_fields)
        fname = agg_fields{i};
        vec = summary.([fname '_vec']);

        % Remove NaN for aggregation
        valid = vec(~isnan(vec));

        if ~isempty(valid)
            summary.([fname '_mean'])   = mean(valid);
            summary.([fname '_std'])    = std(valid);
            summary.([fname '_median']) = median(valid);
            summary.([fname '_min'])    = min(valid);
            summary.([fname '_max'])    = max(valid);
        else
            summary.([fname '_mean'])   = NaN;
            summary.([fname '_std'])    = NaN;
            summary.([fname '_median']) = NaN;
            summary.([fname '_min'])    = NaN;
            summary.([fname '_max'])    = NaN;
        end
    end

    % ================================================================
    % 3. Derived cross-frame metrics
    % ================================================================

    % Contrast stability: CV of contrast_iqr across frames.
    % Low CV = consistent dynamic range (good normalization).
    iqr_vec = summary.contrast_iqr_vec(~isnan(summary.contrast_iqr_vec));
    if ~isempty(iqr_vec) && mean(iqr_vec) > 0
        summary.contrast_stability = std(iqr_vec) / mean(iqr_vec);
    else
        summary.contrast_stability = NaN;
    end

    % Temporal SSIM stability: std of SSIM.
    % Low std = uniformly smooth temporal transitions.
    ssim_vec = summary.temporal_ssim_vec( ...
        ~isnan(summary.temporal_ssim_vec));
    if ~isempty(ssim_vec)
        summary.temporal_ssim_stability = std(ssim_vec);
    else
        summary.temporal_ssim_stability = NaN;
    end

    % SRR stability: std of SRR.
    % For auto-threshold filters, this should be lower.
    srr_vec = summary.srr_vec(~isnan(summary.srr_vec));
    if ~isempty(srr_vec)
        summary.srr_stability = std(srr_vec);
    else
        summary.srr_stability = NaN;
    end

    % Edge sharpness stability: CV of edge sharpness.
    es_vec = summary.edge_sharpness_vec( ...
        ~isnan(summary.edge_sharpness_vec));
    if ~isempty(es_vec) && mean(es_vec) > 0
        summary.edge_sharpness_stability = std(es_vec) / mean(es_vec);
    else
        summary.edge_sharpness_stability = NaN;
    end

    % ================================================================
    % 4. Computational cost
    % ================================================================
    if nargin >= 2 && ~isempty(frame_time)
        ft = frame_time(valid_mask);
        summary.time_mean  = mean(ft);
        summary.time_std   = std(ft);
        summary.time_total = sum(ft);
    else
        summary.time_mean  = NaN;
        summary.time_std   = NaN;
        summary.time_total = NaN;
    end

    % ================================================================
    % 5. Print summary table
    % ================================================================
    fprintf('\n');
    fprintf('============== METRICS SUMMARY ==============\n');
    fprintf('Valid frames:         %d\n', n_valid);
    fprintf('----------------------------------------------\n');
    fprintf('  FILTER SELECTIVITY\n');
    fprintf('    SRR:              %.3f +/- %.3f\n', ...
        summary.srr_mean, summary.srr_std);
    fprintf('    SRR stability:    %.4f\n', ...
        summary.srr_stability);
    fprintf('    Pixel retention:  %.3f +/- %.3f\n', ...
        summary.pixel_retention_mean, ...
        summary.pixel_retention_std);
    fprintf('    Polarity balance: %.3f +/- %.3f\n', ...
        summary.polarity_balance_mean, ...
        summary.polarity_balance_std);
    fprintf('    Clusters:         %.1f +/- %.1f\n', ...
        summary.n_clusters_mean, summary.n_clusters_std);
    fprintf('    Mean cluster sz:  %.1f +/- %.1f px\n', ...
        summary.mean_cluster_size_mean, ...
        summary.mean_cluster_size_std);
    fprintf('----------------------------------------------\n');
    fprintf('  SURFACE QUALITY\n');
    fprintf('    Edge sharpness:   %.4f +/- %.4f\n', ...
        summary.edge_sharpness_mean, ...
        summary.edge_sharpness_std);
    fprintf('    Mean gradient:    %.4f +/- %.4f\n', ...
        summary.mean_gradient_mean, ...
        summary.mean_gradient_std);
    fprintf('    Contrast IQR:     %.4f +/- %.4f\n', ...
        summary.contrast_iqr_mean, ...
        summary.contrast_iqr_std);
    fprintf('    Contrast stab.:   %.4f (CV)\n', ...
        summary.contrast_stability);
    fprintf('    Noise floor:      %.4f +/- %.4f\n', ...
        summary.noise_floor_mean, ...
        summary.noise_floor_std);
    fprintf('    Surface fill:     %.3f +/- %.3f\n', ...
        summary.surface_fill_mean, ...
        summary.surface_fill_std);
    fprintf('----------------------------------------------\n');
    fprintf('  TEMPORAL COHERENCE\n');
    fprintf('    SSIM:             %.4f +/- %.4f\n', ...
        summary.temporal_ssim_mean, ...
        summary.temporal_ssim_std);
    fprintf('    SSIM stability:   %.4f\n', ...
        summary.temporal_ssim_stability);
    fprintf('    MAE:              %.4f +/- %.4f\n', ...
        summary.temporal_mae_mean, ...
        summary.temporal_mae_std);
    fprintf('----------------------------------------------\n');
    fprintf('  COMPUTATIONAL COST\n');
    fprintf('    Time/frame:       %.3f +/- %.3f s\n', ...
        summary.time_mean, summary.time_std);
    fprintf('    Total time:       %.1f s\n', ...
        summary.time_total);
    fprintf('==============================================\n\n');

end