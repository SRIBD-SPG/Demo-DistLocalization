function [coords] = uniformMotion(path, segments, speed, segment_time)
% path: n x 3 matrix representing the path
% segments: the number of time segments to divide the motion into
% speed: the absolute speed of the motion (m/s)
% segment_time: the duration of each segment

if (segments > 0)

    % Calculate the total distance along the path
    distances = vecnorm(diff(path, 1), 2, 2);
    total_distance = sum(distances);

    % Calculate the duration of each time segment
    total_time = segments;
    segment_durations = total_distance / total_time;

    % Calculate the position of each time segment
    segment_positions = [0; cumsum(distances)];
    target_positions = linspace(0, segment_positions(end), segments+1)';

    % Interpolate the path to find the coordinates at each target position
    coords = interp1(segment_positions, path, target_positions, 'linear');

    % Remove the last coordinate, which is duplicated due to the linspace call
    coords = coords(1:end-1,:);
else
    if nargin < 4 
        error('Missing the speed and segment time for calculation.')
    end
    % Earth's radius in meters
    R = 6371000;

    % Convert longitudes and latitudes to radians
    path_rad = path;
    path_rad(:,1:2) = deg2rad(path(:,1:2));

    % Convert path in spherical coordinates to cartesian coordinates
    x = (R + path(:,3)).*cos(path_rad(:,1)).*cos(path_rad(:,2));
    y = (R + path(:,3)).*cos(path_rad(:,1)).*sin(path_rad(:,2));
    z = (R + path(:,3)).*sin(path_rad(:,1));

    % Calculate the total distance along the path using Euclidean distance formula in cartesian coordinates
    distances = sqrt(diff(x).^2 + diff(y).^2 + diff(z).^2);
    total_distance = sum(distances);

    % Calculate the number of segments based on the speed and segment time
    segments = floor(total_distance / (speed * segment_time));

    % Calculate the position of each time segment
    segment_positions = [0; cumsum(distances)];
    target_positions = linspace(0, segment_positions(end), segments+1)';

    % Interpolate the path to find the coordinates at each target position
    cartesian_coords = [x, y, z];
    interpolated_cartesian_coords = interp1(segment_positions, cartesian_coords, target_positions, 'linear');

    % Convert interpolated cartesian coordinates back to spherical coordinates
    r = sqrt(interpolated_cartesian_coords(:,1).^2 + interpolated_cartesian_coords(:,2).^2 + interpolated_cartesian_coords(:,3).^2);
    lat = asin(interpolated_cartesian_coords(:,3) ./ r);
    lon = atan2(interpolated_cartesian_coords(:,2), interpolated_cartesian_coords(:,1));
    alt = r - R;

    % Convert latitudes and longitudes from radians to degrees
    coords_temp = rad2deg([lat, lon]);
    coords = [coords_temp,1000*ones(size(coords_temp,1),1)];

    % Remove the last coordinate, which is duplicated due to the linspace call
    coords = coords(1:end-1,:);
end
end