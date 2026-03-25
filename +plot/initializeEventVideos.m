function [hFigs, hAxs, hImgs, videoWriters] = ...
    initializeEventVideos(cohOut, atsOut, imgSz, videoOutPath)
% INITIALIZEEVENTVIDEOS  Create figure handles and video writers.
%
%   [HFIGS, HAXS, HIMGS, VIDEOWRITERS] = INITIALIZEEVENTVIDEOS(
%   COHOUT, ATSOUT, IMGSZ) initializes off-screen figures and AVI
%   video writers for recording the ATS surface and coherence map
%   outputs during processing.
%
%   Inputs:
%     cohOut - Logical. If true, initialize coherence video output.
%     atsOut - Logical. If true, initialize ATS video output.
%     imgSz  - [1 x 2] Image dimensions [nRows, nCols].
%     videoOutPath - String. Where the output video should be saved.
%
%   Outputs:
%     hFigs        - {1 x 2} Cell array of figure handles.
%     hAxs         - {1 x 2} Cell array of axes handles.
%     hImgs        - {1 x 2} Cell array of image handles.
%     videoWriters - {1 x N} Cell array of VideoWriter objects.
%
%   Notes:
%     - Figures are created with 'Visible', 'off' for headless
%       recording. Call set(hFigs{k}, 'Visible', 'on') to display.
%     - Video paths are set relative to the current directory.
%       Modify the videoFileName variables below if needed.
%     - Video writers are opened upon creation. Call
%       close(videoWriters{k}) when processing is complete.
%
%   See also: VideoWriter

    % ----------------------------------------------------------------
    % 0. Initialize output cells
    % ----------------------------------------------------------------
    hFigs  = {[], []};
    hAxs   = {[], []};
    hImgs  = {[], []};
    videoWriters = {};

    % ----------------------------------------------------------------
    % 1. ATS surface output
    % ----------------------------------------------------------------
    if atsOut
        hFigATS = figure('Visible', 'off', ...
            'Position', [100 100 imgSz]);
        hAxATS  = axes('Parent', hFigATS);
        colormap(hAxATS, 'gray');
        colorbar(hAxATS);
        set(hAxATS, 'FontSize', 16, 'Color', 'white', ...
            'XTick', [], 'YTick', []);
        initial_data = nan(imgSz(2), imgSz(1));
        hImgATS = imagesc(hAxATS, initial_data);
        set(hImgATS, 'AlphaData', ~isnan(initial_data));

        videoFileName = fullfile(videoOutPath);
        videoWriters{1} = VideoWriter(videoFileName);
        videoWriters{1}.FrameRate = 60;
        open(videoWriters{1});

        hFigs{1} = hFigATS;
        hAxs{1}  = hAxATS;
        hImgs{1} = hImgATS;
    end

    % ----------------------------------------------------------------
    % 2. Coherence map output
    % ----------------------------------------------------------------
    if cohOut
        hFigCOH = figure('Visible', 'off', ...
            'Position', [100 100 imgSz]);
        hAxCOH  = axes('Parent', hFigCOH);
        colormap(hAxCOH, 'gray');
        colorbar(hAxCOH);
        set(hAxCOH, 'FontSize', 16, 'Color', 'white');
        initial_data = nan(imgSz(2), imgSz(1));
        hImgCOH = imagesc(hAxCOH, initial_data);
        set(hImgCOH, 'AlphaData', ~isnan(initial_data));

        videoFileName = fullfile(pwd, ...
            'coherence_map_output.avi');
        videoWriters{2} = VideoWriter(videoFileName);
        videoWriters{2}.FrameRate = 60;
        open(videoWriters{2});

        hFigs{2} = hFigCOH;
        hAxs{2}  = hAxCOH;
        hImgs{2} = hImgCOH;
    end

end