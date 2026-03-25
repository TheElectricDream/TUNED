classdef EventBuffer < handle
% EVENTBUFFER  Chunked HDF5 event stream reader.
%
%   BUF = EVENTBUFFER(FILEPATH, T_START, T_END) opens the HDF5 event
%   file at FILEPATH and prepares a buffered reader that streams
%   events in chunks rather than loading the entire recording into
%   memory. Only events in [T_START, T_END] seconds are returned.
%
%   [X, Y, T, P] = BUF.NEXTWINDOW(T_RANGE_N, IMGSZ) returns all
%   events with timestamps up to T_RANGE_N (in shifted time, where
%   t = 0 corresponds to T_START). Events outside [1, imgSz] spatial
%   bounds are filtered. The buffer automatically refills from disk
%   when needed.
%
%   Properties:
%     filepath      - Path to HDF5 file.
%     t_start       - Start of processing window [s] (original time).
%     t_end         - End of processing window [s] (original time).
%     t_total       - Duration of the processing window [s].
%     n_events_file - Total events in the HDF5 file.
%     chunk_size    - Number of events per disk read (default 500000).
%     is_exhausted  - True when all events in [t_start, t_end] have
%                     been consumed.
%
%   Algorithm:
%     1. On construction, read the first timestamp for time
%        normalization and determine the total event count via
%        h5info. No event data is loaded yet.
%     2. On the first call to nextWindow, load the first chunk
%        starting from the file position corresponding to t_start
%        (found via binary search on timestamps).
%     3. On each nextWindow call, scan the buffer for events up to
%        t_range_n. If the buffer is exhausted before t_range_n is
%        reached, load the next chunk from disk and continue.
%     4. Consumed events at the front of the buffer are discarded
%        on each refill to keep memory bounded.
%
%   Notes:
%     - Peak memory is bounded at ~2 * chunk_size * 4 arrays *
%       4 bytes (single) ≈ 16 MB for the default chunk_size.
%     - The HDF5 datasets must be named '/timestamp', '/x', '/y',
%       '/polarity', matching the EVOS / NEXUS convention.
%     - Timestamps in the file are in microseconds (integer).
%       They are converted to seconds and shifted so that
%       t_start maps to t = 0, matching main.m convention.
%
%   Example:
%     buf = process.EventBuffer(fullfile(hdf5Path, fileName), ...
%               t_start_process, t_end_process);
%     while ~buf.is_exhausted
%         [x, y, t, p] = buf.nextWindow(t_range_n, imgSz);
%         % ... process frame ...
%     end
%
%   See also: process.sliceToValidRange, main

    properties
        filepath         char
        t_start          double        % [s] original time axis
        t_end            double        % [s] original time axis
        t_total          double        % [s] duration of window
        n_events_file    int64         % total events in HDF5
        chunk_size       int64 = int64(500000)
        is_exhausted     logical = false
    end

    properties (Access = private)
        t0_file          double        % first timestamp in file [s]
        file_cursor      int64         % next unread row in HDF5 (1-based)
        start_cursor     int64         % file row where t_start begins
        end_cursor       int64         % file row where t_end ends
        buf_x            single
        buf_y            single
        buf_t            single
        buf_p            single
        buf_offset       int64         % local index into buffer
        buf_len          int64         % number of events in buffer
        initialized      logical = false
    end

    methods

        % ============================================================
        % Constructor
        % ============================================================
        function obj = EventBuffer(filepath, t_start, t_end, varargin)
        % EVENTBUFFER  Create a buffered HDF5 event reader.
        %
        %   BUF = EVENTBUFFER(FILEPATH, T_START, T_END)
        %   BUF = EVENTBUFFER(FILEPATH, T_START, T_END, CHUNK_SIZE)

            obj.filepath = filepath;
            obj.t_start  = t_start;
            obj.t_end    = t_end;
            obj.t_total  = t_end - t_start;

            if nargin >= 4
                obj.chunk_size = int64(varargin{1});
            end

            % ---- Read file metadata (no event data loaded) ----
            info = h5info(filepath, '/timestamp');
            obj.n_events_file = int64(info.Dataspace.Size);

            % First timestamp for microsecond-to-second conversion
            t_first = double(h5read(filepath, '/timestamp', 1, 1));
            obj.t0_file = t_first / 1e6;

            % ---- Binary search for start and end cursors ----
            obj.start_cursor = obj.findTimeCursor( ...
                t_start + obj.t0_file);
            obj.end_cursor   = obj.findTimeCursor( ...
                t_end + obj.t0_file);

            obj.file_cursor = obj.start_cursor;

            % Buffer starts empty; first nextWindow triggers fill
            obj.buf_x      = single([]);
            obj.buf_y      = single([]);
            obj.buf_t      = single([]);
            obj.buf_p      = single([]);
            obj.buf_offset = int64(1);
            obj.buf_len    = int64(0);
        end

        % ============================================================
        % nextWindow — replaces sliceToValidRange
        % ============================================================
        function [x_valid, y_valid, t_valid, p_valid] = ...
                nextWindow(obj, t_range_n, imgSz)
        % NEXTWINDOW  Return events up to t_range_n (shifted time).
        %
        %   [X, Y, T, P] = BUF.NEXTWINDOW(T_RANGE_N, IMGSZ)
        %
        %   t_range_n is in shifted seconds (t_start = 0), matching
        %   the convention in main.m.

            % Ensure buffer has data
            if obj.buf_len == 0 || obj.buf_offset > obj.buf_len
                obj.refill();
            end

            if obj.buf_len == 0
                % No events remain
                obj.is_exhausted = true;
                x_valid = []; y_valid = [];
                t_valid = []; p_valid = [];
                return;
            end

            % Scan buffer for events <= t_range_n
            scan_start = obj.buf_offset;
            scan_pos   = obj.buf_offset;

            % We may need to refill mid-scan if the window extends
            % beyond the current buffer
            x_acc = {};
            y_acc = {};
            t_acc = {};
            p_acc = {};

            while true
                % Scan within current buffer
                while scan_pos <= obj.buf_len && ...
                      obj.buf_t(scan_pos) <= t_range_n
                    scan_pos = scan_pos + 1;
                end

                scan_end = scan_pos - 1;

                % Collect events from this buffer segment
                if scan_end >= scan_start
                    x_acc{end+1} = obj.buf_x(scan_start:scan_end); %#ok<AGROW>
                    y_acc{end+1} = obj.buf_y(scan_start:scan_end); %#ok<AGROW>
                    t_acc{end+1} = obj.buf_t(scan_start:scan_end); %#ok<AGROW>
                    p_acc{end+1} = obj.buf_p(scan_start:scan_end); %#ok<AGROW>
                end

                % Did we stop because we hit t_range_n or ran out
                % of buffer?
                if scan_pos <= obj.buf_len
                    % Stopped at t_range_n — done
                    obj.buf_offset = int64(scan_pos);
                    break;
                else
                    % Ran out of buffer — try to refill
                    obj.refill();
                    if obj.buf_len == 0
                        % No more data on disk
                        obj.is_exhausted = true;
                        break;
                    end
                    scan_start = int64(1);
                    scan_pos   = int64(1);
                end
            end

            % Concatenate accumulated segments
            if isempty(x_acc)
                x_valid = []; y_valid = [];
                t_valid = []; p_valid = [];
                return;
            end

            x_all = cat(1, x_acc{:});
            y_all = cat(1, y_acc{:});
            t_all = cat(1, t_acc{:});
            p_all = cat(1, p_acc{:});

            % Spatial bounds filter
            valid_mask = x_all >= 1 & x_all <= imgSz(1) & ...
                         y_all >= 1 & y_all <= imgSz(2);

            x_valid = x_all(valid_mask);
            y_valid = y_all(valid_mask);
            t_valid = t_all(valid_mask);
            p_valid = p_all(valid_mask);
        end

        % ============================================================
        % reset — allow re-processing from the start
        % ============================================================
        function reset(obj)
        % RESET  Rewind the buffer to the start of the processing
        %   window so the frame loop can be re-run without
        %   reconstructing the object.

            obj.file_cursor  = obj.start_cursor;
            obj.buf_x        = single([]);
            obj.buf_y        = single([]);
            obj.buf_t        = single([]);
            obj.buf_p        = single([]);
            obj.buf_offset   = int64(1);
            obj.buf_len      = int64(0);
            obj.is_exhausted = false;
        end

    end

    methods (Access = private)

        % ============================================================
        % refill — load next chunk from HDF5
        % ============================================================
        function refill(obj)
        % REFILL  Read the next chunk_size events from the HDF5 file.

            if obj.file_cursor > obj.end_cursor
                obj.buf_len = int64(0);
                return;
            end

            % How many events remain in the valid range
            remaining = obj.end_cursor - obj.file_cursor + 1;
            n_read    = min(int64(remaining), obj.chunk_size);

            % HDF5 hyperslab read: [start, count]
            fc = double(obj.file_cursor);
            nr = double(n_read);
            
            raw_t = double(h5read(obj.filepath, '/timestamp', fc, nr));
            raw_x = single(h5read(obj.filepath, '/x', fc, nr));
            raw_y = single(h5read(obj.filepath, '/y', fc, nr));
            raw_p = single(h5read(obj.filepath, '/polarity', fc, nr));

            obj.file_cursor = obj.file_cursor + n_read;

            % Convert timestamps: microseconds → seconds, then
            % shift so that t_start = 0
            obj.buf_t = single(raw_t / 1e6 - obj.t0_file ...
                        - obj.t_start);
            obj.buf_x = raw_x + 1;            % 0-indexed → 1-indexed
            obj.buf_y = raw_y + 1;
            obj.buf_p = raw_p;

            obj.buf_offset = int64(1);
            obj.buf_len    = int64(n_read);
        end

        % ============================================================
        % findTimeCursor — binary search for file position
        % ============================================================
        function cursor = findTimeCursor(obj, t_target_sec)
        % FINDTIMECURSOR  Binary search for the first event at or
        %   after t_target_sec in the HDF5 timestamp dataset.

            t_target_us = int64(t_target_sec * 1e6);

            lo = int64(1);
            hi = obj.n_events_file;

            while lo < hi
                mid = lo + idivide(hi - lo, int64(2));
                t_mid = h5read(obj.filepath, '/timestamp', double(mid), 1);
                if t_mid < t_target_us
                    lo = mid + 1;
                else
                    hi = mid;
                end
            end

            cursor = lo;
        end

    end
end