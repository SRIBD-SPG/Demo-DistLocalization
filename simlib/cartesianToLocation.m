function location = cartesianToLocation(x, y, z)
% function to convert Cartesian coordinates to geographic coordinates (latitude, longitude, altitude)
wgs84 = wgs84Ellipsoid('meter');
[lat, lon, alt] = ecef2geodetic(wgs84, x, y, z);
location = [lat, lon, alt]';
end