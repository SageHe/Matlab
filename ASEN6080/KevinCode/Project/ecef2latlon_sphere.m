function LL = ecef2latlon_sphere(Recef,Re)
    lat = asind(Recef(3)/Re);
    lon = acosd(Recef(1)/(Re*cosd(lat)));
    LL = [lat lon];
end