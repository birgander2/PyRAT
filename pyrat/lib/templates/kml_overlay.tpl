<%
  import numpy as np
  import os
%>
<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2">
<Document>
<GroundOverlay>
  <name>${name}</name>
      <Icon>
        <href>${os.path.basename(png_file)}</href>
     </Icon>
     <gx:LatLonQuad>
     <coordinates>
        ${' '.join(['%.12f,%.12f'%(np.rad2deg(c[0]),np.rad2deg(c[1])) for c in corners])}
     </coordinates>
     </gx:LatLonQuad>
</GroundOverlay>
</Document>
</kml>
