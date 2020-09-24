function gps = GPS_newparser(file_name,outputfile_name)
% GPS_parser.m: This code takes in the file name for raw NEMA output from the 
% adafruit ulimate GPS sensor, looks for the GGA sentence and takes the 
% lat, long, alt and converts it to degrees.
%
% Tyler Clayton
% 4/15/2016
%
% INPUTS:   file_name -         sting of text file name
%           outputfile_name -   name of output text file
%
% OUTPUTS:  gps       -  structure containing:
%                          gps.lat, gps.long - latitude/longitude [degrees]
%                          gps.alt           - altitude [meters]
%                          gps.totpts        - total points taken
%                          gps.validSamp        - num of valid gps points
%                          gps.samplesB4fix  - # of void samples recorded
%                                              before initial fix aquired
%                          gps.voidSamp      - # of void samples after
%                                              after initial fix
%%%

%% Parse the GPS data 

% open text file with gps data  
fid_in = fopen(file_name,'r');
% create test file for output
fid_out = fopen(outputfile_name,'w');
%add header to file
% if(strcmp(output_type,'WP')) % for waypoint output file, not a continuous track
%     header_str = 'latitude,longitude,altitude\n';
% else %default to track
    header_str = 'type,latitude,longitude,altitude\n';
% end
fprintf(fid_out,header_str);


f = 0; % counter for # of fix
k = 0; % counter for # of samples

while 1
    gps_unformat = fgetl(fid_in);
    %check if the first string read in is empty. if so take the next one
    if(k==0 && isempty(gps_unformat))
        gps_unformat = fgetl(fid_in);
    end
    
    if ~ischar(gps_unformat), break, end
    % parse the data, look for correct NEMA sentence -> $GPGGA -> lat/long/alt
    if strcmp(gps_unformat(1:6), '$GPGGA') %gps_unformat(1:6) == '$GPGGA'
        %if greater string>38 chars has gps data, otherwise no fix was found
        if(~contains(gps_unformat,',,,,,'))
            if(f==0)
                gps.samplesB4fix = k; % mark the initial number of samples recorded before a fix was aquired
            end
            f = f+1; %increment fix counter
            % need to change lat/long data from ddmm.mmmm to degrees only
            gps.t.hr(f) = str2double(gps_unformat(8:9));
            gps.t.min(f) = str2double(gps_unformat(10:11));
            gps.t.s(f) = str2double(gps_unformat(12:17));
            if(f==1)
               t_init = gps.t.hr(f)*3600 + gps.t.min(f)*60 + gps.t.s(f);
               dNum = 0;
            end
            gps.time(f) = (gps.t.hr(f)*3600 + gps.t.min(f)*60 + gps.t.s(f)) - t_init;
            if(gps.time(f)<0)
               gps.time(f) = gps.time(f) + 3600*24*dNum;
               dNum = dNum+1;
            end
            gps_lat_min = str2double(gps_unformat(21:27));
            gps_lat_deg = str2double(gps_unformat(19:20));
            gps.lat(f) = gps_lat_deg + gps_lat_min/60;
            if(strcmp(gps_unformat(29),'S'))
                gps.lat(f)=-1*gps.lat(f);
            end
            
            gps_long_min = str2double(gps_unformat(34:40));
            gps_long_deg = str2double(gps_unformat(31:33));
            gps.long(f) = gps_long_deg + gps_long_min/60;
            if(strcmp(gps_unformat(42),'W'))
                gps.long(f)=-1*gps.long(f);
            end

            gps.alt(f) = str2double(gps_unformat(54:59));
            
            %%% WRITE TO TXT FILE FOR GPS VISUALIZER
%             if(strcmp(output_type,'WP')) % for waypoint output file, not a continuous track
%                 fprintf(fid_out, '%2.4f,%3.4f,%6.2f\n',gps.lat(f),gps.long(f),gps.alt(f));
%             else %default to track
                fprintf(fid_out, 'T,%2.4f,%3.4f,%6.2f\n',gps.lat(f),gps.long(f),gps.alt(f));
%             end
            
            
            
        end
        k=k+1; %increment sample counter
    end
end
fclose(fid_in);
fclose(fid_out);
gps.totpts = k;
gps.validSamp = f;
gps.voidSamp = gps.totpts-(gps.validSamp+gps.samplesB4fix);
end

