function position = locationToCartesian(lat, lon, alt)
% function to convert geographic coordinates (latitude, longitude, altitude) to Cartesian coordinates
wgs84 = wgs84Ellipsoid('meter');
[x, y, z] = geodetic2ecef(wgs84, lat, lon, alt);
position = [x, y, z]';

end