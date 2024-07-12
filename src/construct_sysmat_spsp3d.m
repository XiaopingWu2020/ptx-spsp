function sysmat = construct_sysmat_spsp3d (grad, b1maps, mask, fox, b0map, freqs, ...
                                           wts, tw,dt,poffset,rftype)
% CONSTRUCT_SYSMAT_SPSP3D Construct the system matrix for designing 2D spatial and 1D spectral pTX RF.
%
% Usage: sysmat = construct_sysmat_spsp3d (grad, b1maps, mask, fox, b0map, freqs, wts, 
% tw,dt,poffset,rftype)
%
% Returns
% -------
% sysmat: 
%
% Expects
% -------
% grad: 
% b1maps: 
% mask: 
% fox: 
% b0map: [] for zero.
% freqs: a vector of chemical shift related freq offsets in hz that are of
% interest.
% wts: a vector of weights used to control the relative importance of freqs.
% defaults to a unity vector, i.e., ones(size(freqs)).
% tw: time window
% dt: dwell time
% poffset: spatial offset.
% 
% rftype: specifies what linear type of rf pulse will be designed.0 for small
% tip angle (sta) pulse design and 1 for linear class large tip angle (lclta)
% pulse design. defaults to 1.
%
%
% See also: construct_system_matrix construct_sysmat_lclta
% construct_sysmat_spspkT
%
%
% Copyright (C) 2010 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Fri Apr 16 15:43:59 2010
%

if nargin < 11
  rftype = 1;
end
if nargin<10
    poffset=[0 0 0];
end
if nargin<9
    dt=10e-6;
end
if nargin<8
    tw=[];
end

if nargin < 5 || isempty(b0map)
  b0map = zeros(size(mask));
end

if nargin< 7|| isempty(wts)
  wts = ones(size(freqs));
end

%
nfreqs = length(freqs);

nchs = size(b1maps,4);
nt = size(grad,2);
nspts = length(mask(mask));

if isempty(tw)
sysmat = complex(zeros(nfreqs*nspts,nchs*nt));else
sysmat = complex(zeros(nfreqs*nspts,nchs*length(tw(tw))));
end


switch rftype
  case 0
  for ifreq= 1:nfreqs,
    iIdx0 = (ifreq-1)*nspts + 1;
    isysmat = construct_system_matrix(grad, b1maps, mask, ...
                                                      fox, b0map+hz2tesla(freqs(ifreq)),...
                                                      tw,[],dt,poffset);
      sysmat(iIdx0:iIdx0+nspts-1,:)= wts(ifreq).* isysmat;                                            
  end
 
 case 1
  for ifreq= 1:nfreqs,
    iIdx0 = (ifreq-1)*nspts + 1;
    sysmat(iIdx0:iIdx0+nspts-1,:) = wts(ifreq).* construct_sysmat_lclta(grad, b1maps, mask, ...
                                                      fox, b0map+ ...
                                                      hz2tesla(freqs(ifreq)),tw,dt);
  end
 otherwise
end


disp('-> System matrix ready for SPSP pulse design!');


