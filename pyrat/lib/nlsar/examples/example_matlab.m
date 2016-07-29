%%
%% example_matlab.m: Example of a MATLAB program using NL-SAR Toolbox
%%
%% This file is part of NL-SAR Toolbox version 0.6.
%%
%% Copyright Charles-Alban Deledalle (2013)
%% Email charles-alban.deledalle@math.u-bordeaux1.fr
%%
%% This software is a computer program whose purpose is to provide a
%% suite of tools to manipulate SAR images.
%%
%% This software is governed by the CeCILL license under French law and
%% abiding by the rules of distribution of free software. You can use,
%% modify and/ or redistribute the software under the terms of the CeCILL
%% license as circulated by CEA, CNRS and INRIA at the following URL
%% "http://www.cecill.info".
%%
%% As a counterpart to the access to the source code and rights to copy,
%% modify and redistribute granted by the license, users are provided only
%% with a limited warranty and the software's author, the holder of the
%% economic rights, and the successive licensors have only limited
%% liability.
%%
%% In this respect, the user's attention is drawn to the risks associated
%% with loading, using, modifying and/or developing or reproducing the
%% software by the user in light of its specific status of free software,
%% that may mean that it is complicated to manipulate, and that also
%% therefore means that it is reserved for developers and experienced
%% professionals having in-depth computer knowledge. Users are therefore
%% encouraged to load and test the software's suitability as regards their
%% requirements in conditions enabling the security of their systems and/or
%% data to be ensured and, more generally, to use and operate it in the
%% same conditions as regards security.
%%
%% The fact that you are presently reading this means that you have had
%% knowledge of the CeCILL license and that you accept its terms.
%%
%%
%% Started on  Fri Aug 23 18:14:48 2013 Charles-Alban Deledalle
%% Last update Sat Aug 24 13:32:38 2013 Charles-Alban Deledalle
%%

close all;
clear all;

try
  [M, N, D] = sarinfo('example.rat')

  sarimage = sarread('example.rat');

  sarimage_box   = sarboxcar(sarimage);
  sarimage_disk  = sardiskcar(sarimage);
  sarimage_gauss = sargausscar(sarimage);
  [sarimage_nl, eqlook_nl] = sarnlsar(sarimage, 3);

  figure;
  subplot(2, 3, 1);
  sarshow(sarimage);
  subplot(2, 3, 2);
  sarshow(sarimage_nl);
  subplot(2, 3, 3);
  imshow(eqlook_nl'); colormap(gray); axis image;
  caxis([0, 100]);
  axis off;
  subplot(2, 3, 4);
  sarshow(sarimage_box);
  subplot(2, 3, 5);
  sarshow(sarimage_disk);
  subplot(2, 3, 6);
  sarshow(sarimage_gauss);
  linkaxes;

catch error
  disp(error.message);
  disp('Please run from command line:');
  disp('       sarmire example.rat 256 256 3 3');

end