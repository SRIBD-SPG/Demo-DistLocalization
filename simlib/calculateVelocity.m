function avg_velocity = calculateVelocity(start_point, end_point, time)

% load assistant function
if exist('locationToCartesian', 'file') ~= 2
    addpath(fullfile('..', '..', '..', 'ShareFolder', 'Commons'));
end

avg_velocity = [];

% 计算距离差
for k = 1:size(start_point, 1)

    delta_distance = locationToCartesian(end_point(k, 1), end_point(k, 2), end_point(k, 3)) ...
        - locationToCartesian(start_point(k, 1), start_point(k, 2), start_point(k, 3));

    % 计算平均速度
    avg_velocity = [avg_velocity; delta_distance' / time];
end