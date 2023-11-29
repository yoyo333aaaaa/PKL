
%------------------------------------------------------------------------
%
% ggmplus2013_v4
%
%------------------------------------------------------------------------
% Purpose:  Seamless extraction of GGMplus gravity field functionals
%           from GGMplus binary files
%------------------------------------------------------------------------
%
% Input:
%    functional:  GGMplus functional to extract.
%                 values: 'geoid', 'gravity','acceleration' 'xi', 'eta'
%
%    minlon1/maxlon1/minlat1/maxlat1: geodetic coordinates of 
%                                     the selection area
%    
%    facX,facY:   resolution scaling factors
%                 facX,facY = 1 : returns data at 7.2 arc-sec resolution
%                                 no interpolation.
%                 facX,facY > 1 : upsamples the data to higher resolution
%                                 than 7.2 arc-sec, interpolation required
%                 facX,facY < 1 : downsamples the data to lower resolution
%                                 than 7.2 arc-sec, interpolation required
%    path_basis:       Name of path where the GGMplus data folders are located
%    
%    nodataval:   Desired value to be assigned outside the data area
%                 (e.g, offshore or North of 60 deg/ South of -60 deg
%                 latitude).
% 
%    method:      interpolation method applied when facX,facY are not   
%                 equal to 1. 'cubic' is recommended. Parameter is
%                 ignored when facX,facY = 1.
%
%------------------------------------------------------------------------
%  Output:       
%                X  meshgrid-compatible matrix with longitudes
%                Y  meshgrid-compatible matrix with latitudes
%                Z  matrix with extracted GGMplus functional
%
%------------------------------------------------------------------------
%
% Example call:
%
% [X,Y,Z] = ggmplus2013_v4('gravity', 110.00,110.75,-8.00,-7.50,1,1,'D:\PKL\Gravity',0,'linear')
%
%------------------------------------------------------------------------
%
% Notes:
%
% -when facX=facY=1, data is extracted without interpolation.
%  This means that X and Y do not necessarily match 
%  exactly the minlon1/maxlon1/minlat1/maxlat1 selection area
%
% -when there is interpolation, data is interpolated to a line-registrated
%  grid matching exactly  minlon1,maxlon1,minlat,maxlat1
%
%------------------------------------------------------------------------
% Christian Hirt
% WA Centre for Geodesy
% Curtin University
% Perth
% last edited 2013-05-30
%------------------------------------------------------------------------

function [X,Y,Z] = ggmplus2013_v4(functional, minlon1,maxlon1,minlat1,maxlat1,...
                           facX,facY,path_basis,nodataval,method)

                      
curr_dir = pwd;

%-------------------------------------------------------------------------
% 1 - assign constant parameters (res, tilesize, infofile)
%-------------------------------------------------------------------------

% some constant that apply for all GGMplus grids -------------------------
 
machineform = 'ieee-be';
tilesize    = 5;
resolution  = 7.2/3600;
tileelems   = 2500;          % elements in lat and lon direction


% constants that are different for the functionals ------------------------
switch functional
    case {'geoid'}    
       conv_factor = 1000;  % mm to m
       suffix      = '.ha';
       path_tiles  = fullfile(path_basis, 'geoid'  ); 
       storagetype = 'int32';
       init_nodata   = -2^31;

    case {'gravity','grav_disturbance'}   
       conv_factor = 10;  % 0.1 mGal to mGal
       suffix      = '.dg';
       path_tiles  = fullfile(path_basis, 'dg'  ); 
       storagetype = 'int16';
       init_nodata   = -2^15;
 
    case {'acceleration','grav_acceleration'}   
       conv_factor = 10;  % 0.1 mGal to mGal
       suffix      = '.ga';
       path_tiles  = fullfile(path_basis, 'ga'  ); 
       storagetype = 'int32';
       init_nodata   = -2^31;
  
    case {'xi'}    
       conv_factor = 10;  % 0.1 arc-sec to arc-sec
       suffix      = '.xi';
       path_tiles  = fullfile(path_basis, 'xi'  ); 
       storagetype = 'int16';
       init_nodata   = -2^15;
  
    case {'eta'} 
       conv_factor = 10;  % 0.1 arc-sec to arc-sec
       suffix      = '.eta';
       path_tiles  = fullfile(path_basis, 'eta'  );       
       storagetype = 'int16';
       init_nodata   = -2^15;
end


% * fac 1 means no interpolation carried out, data is passed through.
% * fac > 1 or < 1 means interpolation required and selection area 
%   extended by a small strip of 0.1 deg width

if(facX==1&facY==1)  % no interpolation
    interpflag = 0;
    tolarea    = 0;
else
    interpflag = 1;
    tolarea    = 0.1;   % extension in deg to either side 
  
end


%------------------------------------------------------------------------
% 2 extend selection area to match the dimensions of dem tiles
%------------------------------------------------------------------------

% minlon1,maxlon1 etc are the original boundaries passed to the function
% minlon2,maxlon2 etc are those extended by areatol. If facs ==1, the
% values are the same

minlon2=minlon1-tolarea;
maxlon2=maxlon1+tolarea;
minlat2=minlat1-tolarea;
maxlat2=maxlat1+tolarea;

% min coordinates of most eastern, and western etc tiles to be accessed
Sall_minlon = floor (minlon2/tilesize)*tilesize;
Sall_maxlon = floor (maxlon2/tilesize)*tilesize;
Sall_minlat = floor (minlat2/tilesize)*tilesize;
Sall_maxlat = floor (maxlat2/tilesize)*tilesize;

Sall_lonvec =Sall_minlon:tilesize:Sall_maxlon;
Sall_latvec =Sall_minlat:tilesize:Sall_maxlat;
sn=length(Sall_lonvec);
sm=length(Sall_latvec);

[Sall_lonmat Sall_latmat]  = meshgrid(Sall_lonvec,Sall_latvec);

Slonmat= reshape(Sall_lonmat,sn*sm,1);
Slatmat= reshape(Sall_latmat,sn*sm,1);
ntiles = length(Slonmat);
clear Sall_*

%------------------------------------------------------------------------
% 3 get filenames for the tiles that will be accessed
%------------------------------------------------------------------------
no = 0;
for min_lon = -180-tilesize:tilesize:180
    for min_lat = -90:tilesize:90-tilesize

    no = no +1;
    minlonf = min_lon;
    minlatf  = min_lat;
    if(min_lon<-180)
      minlonf = minlonf+360;
    end
    if(min_lon>180-tilesize)
      minlonf = minlonf-360;
    end
      
    temp1 = getfilename2012(minlatf,minlonf);
    outfilename = strcat(temp1,suffix);

    fchar(no,1) = cellstr(outfilename);
    flon(no) = min_lon;
    flat(no) = min_lat;

  end
end

flat = flat';
flon = flon';
fchar = char(fchar);

%------------------------------------------------------------------------
% 4 intersect
%------------------------------------------------------------------------
ix= [];
for i = 1:ntiles
    ixt = find(flat==Slatmat(i)&flon==Slonmat(i));
    ix = [ix ixt];
end
% ix is now a vector containing indices to all source tiles

% these are the tiles that we have to access
slatmin =flat(ix);              % vectors (ntiles, 1)
slonmin =flon(ix);
slatmax =slatmin+tilesize;
slonmax =slonmin+tilesize;

disp(' tiles to be accessed')
tilenames = char(fchar(ix,:));  % vectors (ntiles, 1)

ntilesfound = length(ix);

clear Slatmat Slonmat fchar flat flon

%------------------------------------------------------------------------
% 5 initialise target matrix
%------------------------------------------------------------------------

% use integer coordinates to 0.01 sec cell-centred
fac=3600*100;
r = int32(resolution*fac);  % resolution 
r2 =int32(r/2);             % half-resolution

% this determines the definite geometry of the target matrix
Amin_lat = int32(floor(minlat2*fac/r)*r+r2);
Amax_lat = int32(floor(maxlat2*fac/r)*r-r2);
Amin_lon = int32(floor(minlon2*fac/r)*r+r2);
Amax_lon = int32(floor(maxlon2*fac/r)*r-r2);
% corner coordinates and vectors and r in 0.01 arc sec and exactly match
% the source data

Avec_lat = Amin_lat:r:Amax_lat;
Avec_lon = Amin_lon:r:Amax_lon;
An1 = length(Avec_lat);
Am1 = length(Avec_lon);
  
% Initialize target tile
clear A
A=double(ones(An1,Am1)*NaN);

%------------------------------------------------------------------------
% 6 fill the matrix with data from binary files
%------------------------------------------------------------------------
for i = 1:ntilesfound

     %  step A ----------------------------------------------------------
     %  build vectors for the source matrix
     Smin_lat =  int32(slatmin(i)*fac+r2);
     Smax_lat =  int32(slatmax(i)*fac-r2);
 
     Smin_lon =  int32(slonmin(i)*fac+r2);
     Smax_lon =  int32(slonmax(i)*fac-r2);
    
     Svec_lat =  Smin_lat:r:Smax_lat;
     Svec_lon =  Smin_lon:r:Smax_lon;
     clear Smin_lat Smax_lat Smin_lon Smax_lon
    
     % step B intersect and return common lat/lon values ----------------
     a_lat = intersect (Svec_lat, Avec_lat);
     a_lon = intersect (Svec_lon, Avec_lon);
    
     n_elems_lat = length(a_lat);
     n_elems_lon = length(a_lon);
     if(n_elems_lat&n_elems_lon) % both have be at least one
     
     % step C find the four source and four target indices --------------
     % Destination indices lat
     Aix_lat_min = find(Avec_lat==a_lat(1));
     Aix_lat_max = find(Avec_lat==a_lat(end));
     % Source indices lat
     Six_lat_min = find(Svec_lat==a_lat(1));
     Six_lat_max = find(Svec_lat==a_lat(end));
       
     % Destination indices lat
     Aix_lon_min = find(Avec_lon==a_lon(1));
     Aix_lon_max = find(Avec_lon==a_lon(end));
     % Source indices lat
     Six_lon_min = find(Svec_lon==a_lon(1));
     Six_lon_max = find(Svec_lon==a_lon(end));
        
     % step D open file only if exists -----------------------------------
     te = tileelems;      
     tilename_in = tilenames(i,:);
     %a = strcat('read file',{'   '},tilename_in);
     disp('read file')
     disp(tilename_in)
     cd(path_tiles)   
     if(exist(tilename_in))
            
        fid = fopen(tilename_in,'rb');
        S = fread(fid,te*te,storagetype,machineform);
        fclose(fid);
        disp('read complete')
        S    = reshape(S,te,te);
                   
        %%CH 20130527 temporary use of NaN over no-data areas to avoid
        % interpolation edge effects 
        S(find(S<init_nodata+1))=NaN;          
   
        % copy data block
        A(Aix_lat_min:Aix_lat_max,Aix_lon_min:Aix_lon_max)= ... 
        S(Six_lat_min:Six_lat_max,Six_lon_min:Six_lon_max);

   
        clear  Svec_lon  Svec_lat Aix* Six*  a_lat a_lon S
       
     else
         disp(['file ',tilename_in,' does not exist'])       
     end
     else
         disp(' this tile does not intersect with target area')
     end
end

length(find(isnan(A)))


%------------------------------------------------------------------------
% 7 build the X and Y matrices
%------------------------------------------------------------------------
% convert the basis vectors from integer and 0.01 sec to double and deg
Avec_lon = double(Avec_lon)/fac;
Avec_lat = double(Avec_lat)/fac;

[AX AY]  = meshgrid(Avec_lon,Avec_lat);


%------------------------------------------------------------------------
% 8 interpolate only if fac not 1
%------------------------------------------------------------------------
    
cd(curr_dir)

if(interpflag == 1)
   disp(strcat(' interpolate',{' '},method))
   
   % build target grid, exact match with input minmax boundaries ---------
   res_targetX = resolution*facX; 
   res_targetY = resolution*facY;
   
   % use exact area (minlon1, maxlon1 etc) for interpolation
   B_veclon = minlon1:res_targetX:maxlon1;
   B_veclat = minlat1:res_targetY:maxlat1;
  
   [BX BY]  = meshgrid(B_veclon,B_veclat);

   % carry out the interpolation in blocks -------------------------------   
   
   [n2 m2] = size(BX);
   anz_bloecke = 5;
   blocksize2 = floor(n2/anz_bloecke);
   min2 = 1;
   max2 = min2+blocksize2;

   for i = 1:anz_bloecke
      B(min2:max2,1:m2) = interp2(AX,AY,A,BX(min2:max2,1:m2),BY(min2:max2,1:m2),method);
      min2 = min2+blocksize2+1;
      max2 = min2+blocksize2;
   
      if (i==anz_bloecke-1)
        max2 = n2;
      end    
   end
   
   clear AX AY A
   X = BX;
   clear BX
   Y = BY; 
   clear BY
   Z = B;
   clear B
   disp(' interpolation finished')
 
else 
   disp(' pass extracted data directly.')
 
   X = AX;
   clear AX
   Y = AY; 
   clear AY
   Z = A;
   clear A
       
end

%------------------------------------------------------------------------
% 9 apply scale factor
%------------------------------------------------------------------------
Z = Z/conv_factor; % this scales the GGMplus matrix to the basic unit


%------------------------------------------------------------------------
% 10 find all nodata values and assign user-defined nodataval
%------------------------------------------------------------------------
Z(isnan(Z))=nodataval;


cd(curr_dir)
clear A* B* C* D* E* F* G* H* I* L* M* N* P* T*
clear a* b* c* d* e* f* g* h* i* l* m* n* p* t* s* r*




%------------------------------------------------------------------------
% Function getfilename2012
%------------------------------------------------------------------------
function name = getfilename2012(lat, lon)

d_int = lat;
c_int = lon;


if (d_int>=0)
      temp1 = 'N';
  else
      temp1 = 'S';
  end
  if (abs(d_int)<10)
    temp1 = strcat(temp1,'0');
  end
 
  c = sprintf('%i',c_int);
  d = sprintf('%i',abs(d_int));

  temp3 = strcat(temp1,d);
  
  if (c_int>=0)
    temp4 = strcat(temp3,'E');
  else
    temp4 = strcat(temp3,'W');
  end
  
  if (abs(c_int)<100)
    temp4 = strcat(temp4,'0');
  end
  if (abs(c_int)<10)
    temp4 = strcat(temp4,'0');
  end
  c = sprintf('%i',abs(c_int));
  name = strcat(temp4,c);

end

end

%------------------------------------------------------------------------
% End
%------------------------------------------------------------------------

