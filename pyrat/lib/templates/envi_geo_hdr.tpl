<%
  import os
  import numpy as np
%>ENVI
description         = {Geocoded F-SAR data: ${os.path.basename(file)}}
samples             = ${hdr.Rat.idl_shape[0]}
lines               = ${hdr.Rat.idl_shape[1]}
bands               = 1
map info            = {UTM, 1, ${hdr.Rat.idl_shape[1]}, ${hdr.Geo.min_east}, ${hdr.Geo.min_north}, ${hdr.Geo.ps_east}, ${hdr.Geo.ps_north}, ${np.abs(hdr.Geo.zone)}, ${'North' if hdr.Geo.zone >= 0 else 'South'},WGS-84,units=Meters}
data type           = ${hdr.Rat.var}
interleave          = BIP
sensor type         = SAR
byte order          = 0
header offset       = 1000
file type           = ENVI Standard
spheroid_name       = WGS-84
georeferenced       = UTM
projection_zone     = ${abs(hdr.Geo.zone)}
pixel_spacing_east  = ${hdr.Geo.ps_east}
pixel_spacing_north = ${hdr.Geo.ps_north}
min_easting         = ${hdr.Geo.min_east}
min_northing        = ${hdr.Geo.min_north}
max_easting         = ${hdr.Geo.min_east+(hdr.Rat.idl_shape[0]-1)*hdr.Geo.ps_east}
max_northing        = ${hdr.Geo.min_north+(hdr.Rat.idl_shape[1]-1)*hdr.Geo.ps_north}
