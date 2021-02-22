function ecef = latlon2ecef(latlon,angle_type,Re)
    % for a spherical Earth

    lat = latlon(1);
    lon = latlon(2);
    
    if angle_type == 'degrees'
        x = cosd(lat)*cosd(lon)*Re;
        y = cosd(lat)*sind(lon)*Re;
        z = sind(lat)*Re;
        ecef = [x y z];
    elseif angle_type == 'radians'
        x = cos(lat)*cos(lon)*Re;
        y = cos(lat)*sin(lon)*Re;
        z = sin(lat)*Re;
        ecef = [x y z];
    end
end