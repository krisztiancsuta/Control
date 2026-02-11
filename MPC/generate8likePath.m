function r = generate8likePath(t, max_yaw_rate)
    % Function to generate yaw rate reference for a simple figure-8 path
    % t: array of time points
    % max_yaw_rate: the maximum yaw rate for the circular sections
    
    r = zeros(size(t)); 
    
    % Split the time vector into sections (straight, curve, straight, curve)
    num_segments = 4;  % Two straight lines and two curves
    segment_length = floor(length(t) / num_segments);  % Length of each segment
    
    % First straight section (yaw rate = 0)
    start_idx = 1;
    end_idx = segment_length;
    r(start_idx:end_idx) = 0;  % Zero yaw rate for straight path
    
    % First curve section (positive yaw rate)
    start_idx = end_idx + 1;
    end_idx = 2 * segment_length;
    r(start_idx:end_idx) = max_yaw_rate;  % Constant positive yaw rate
    
    % Second straight section (yaw rate = 0)
    start_idx = end_idx + 1;
    end_idx = 3 * segment_length;
    r(start_idx:end_idx) = 0;  % Zero yaw rate for straight path
    
    % Second curve section (negative yaw rate)
    start_idx = end_idx + 1;
    end_idx = length(t);
    r(start_idx:end_idx) = -max_yaw_rate;  % Constant negative yaw rate
end
