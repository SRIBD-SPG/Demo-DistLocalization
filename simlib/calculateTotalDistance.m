function total_distance = calculateTotalDistance(path)
    % 地球平均半径（单位：米）
    R = 6371e3;
    
    % 初始化总距离
    total_distance = 0;
    
    % 将经纬度从度转换为弧度
    path_rad = deg2rad(path);
    
    % 计算每一段的距离，并累加
    for i = 1:size(path, 1)-1
        % 计算经纬度的差值
        delta_lat = path_rad(i+1, 1) - path_rad(i, 1);
        delta_lon = path_rad(i+1, 2) - path_rad(i, 2);
        
        % 使用 Haversine 公式计算地面距离
        a = sin(delta_lat/2)^2 + cos(path_rad(i, 1)) * cos(path_rad(i+1, 1)) * sin(delta_lon/2)^2;
        c = 2 * atan2(sqrt(a), sqrt(1-a));
        d = R * c;
        
        % 计算高度差
        delta_h = path(i+1, 3) - path(i, 3);
        
        % 使用勾股定理计算总距离
        total_distance = total_distance + sqrt(d^2 + delta_h^2);
    end
end