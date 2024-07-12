%% This is a demo script showing how to design ptx spatial spectral pulse (2D space + 1D spectral) that can trace dB0. 
%
% Created and edited by Xiaoping, 7/12/2024

close all; clearvars
%% prep
load calibrationMS.mat    % load in ptx calibration 
nslices= size(b0mapMS,3);
soi= round(0.5*nslices);
b1maps= 1e-6*conj(rfmapMS(:,:,soi,:)); %b1mapsMSn(:,:,soi,:);
b0map= 1e-6*(-1).* b0mapMS(:,:,soi); % due to a left hand system. 
mask= maskMS(:,:,soi);img= imgb0(:,:,soi);
%%

fox =  1e-3*(210*[1 1]); %1e-3*(192*[1 1]); %1e-3*[256*[1 0.688]];% 1e-3*192*[1 1]; % in m

nchs = size(b1maps,4);
poffset= [-30 0 0];%[-20 0 0]; % mm, minus for shift to A
dt=10e-6;

%figure, myimagesc(sum(abs(b1maps),4), mask), caxis auto
%figure, myimagesc(img, mask), caxis auto
%h= imfreehand;
%targ0= h.createMask;
%targ= imfilter(double(targ0),fspecial('gauss'));
%figure, myimagesc(targ,mask), caxis([0 1]) %,daspect([1 0.688 1])
%save targ targ targ0

%load targ

%%
myrf=load('rfpat_subpulse.mat');
gradbody= [[myrf.gy_T_m];[myrf.gx_T_m]];
grad0= gradbody;
tw0= logical(myrf.twin);

% % usage: [k,g,s,time] = design_toptSpiral(Nitlv, isRotationallyVariant, res, fov, radius, safetyMargin)
% R= 4;
% res= max(fox)./max(size(mask));
% [~,g] = design_toptSpiral(R, false, 4.*res, [max(fox) max(fox)], [0 1], 0.9);
% grad0= 1e-3.* fliplr(g(:,1:2).');
%%
% mygz=load('37segments');
% gzcomp= 2*mygz.gzcomp;
% hsw0= mygz.rf;

%%
% design spectral select
fatbnd = (-1300:125:-800); %+ 100;
waterbnd = -250:125:250;
freqs = [fatbnd waterbnd];
wsb = 8;
wpb = 1;
% water imaging
spect = [zeros(size(fatbnd)) ones(size(waterbnd))];
wts = [wsb*ones(size(fatbnd)) wpb*ones(size(waterbnd))];
figure, plot(freqs,spect,'x-')
xlabel('Chemical shift (Hz)')
ylabel('Desired excitation (a.u.)')

%% 
nSubRfs= 12;
%timewin= false(1,nSubRfs*length(tw0));
tw= repmat(tw0,[1 nSubRfs]);
grad= repmat(grad0,[1 nSubRfs]);

sysmat = construct_sysmat_spsp3d(grad,b1maps,mask,fox,b0map, -freqs,wts,tw,dt,poffset,0);
%%
fa= 10;
targ1= targ; 
targvect = fa.* construct_targvect_spsp (targ1, mask, spect);

rf = calc_rf_cgls(sysmat,targvect,nchs,tw);
%%%
%rfobj= rfPulse(rf);
%rfobj.TimeStep=dt;
%figure, rfobj.plot_amp
%figure, plot(sum(abs(rf),1))

%grdobj= gradPulse(grad);
%grdobj.TimeStep= dt;
%figure, grdobj.plot

%%
%figure, 
%subplot(1,3,1), plot(freqs,spect,'x-')
%xlabel('Chemical shift (Hz)')
%ylabel('Desired excitation (a.u.)')
%title('Excitation target')
%subplot(1,3,2), rfobj.plot_amp
%subplot(1,3,3), grdobj.plot

%% write ini file
[Nc, Nt]= size(rf);
myrf= reshape(rf', [Nc*Nt, 1]);
grad(3,:)= 0;
%%%
mygrad= grad;
orient= 'sagA2P'; %'transA2P';
switch orient
    case 'transA2P'
        mygrad(1,:)= grad(2,:);
        mygrad(2,:)= -grad(1,:);
    case 'sagA2P'
        mygrad(1,:)= grad(3,:); 
        mygrad(2,:)= -grad(1,:);
        mygrad(3,:)= -grad(2,:);
        
    otherwise
end

mygrad= 1e3* mygrad.';

opt.NOMFLIPANGLE= fa;
opt.FACTOROVERSAMPLE= 1;
opt.RFPULSE_COMMENT= 'spsp';opt.VERBOSE            = false;
save_pTXRFPulse_toINI( mygrad, myrf, [], opt);

%%
errs= zeros(size(freqs));
ccs=errs;
clear mxypatptx2d
%grad1= switch_grad_polarity(grad,[1,-1,1]);
for idx=1:length(freqs)
imxypatptx2d = run_bloch_sim ((rf), grad(1:2,:),(b1maps),mask,fox,b0map+hz2tesla(-freqs(idx)),...
    0,[],dt,poffset);
% imxypatptx2d = run_bloch_sim ((rf), grad(1:2,:),(b1maps),mask,fox,b0map,...
%     0,[],dt,poffset);

% imxypatptx2d = run_bloch_sim ((rf),grad,(b1maps),mask,fox,...
%     zeros(size(b0map))+hz2tesla(freqs(idx)),...
%     0,[],dt,poffset);

errs(idx)= norm(targ(mask)- abs(imxypatptx2d(mask)));
ccs(idx)= corr(targ(mask), abs(imxypatptx2d(mask)));
mxypatptx2d(:,:,idx)= imxypatptx2d;
end

%figure, myimagesc(asind(abs(mxypatptx2d(:,:,end))),mask), colormap jet,axis square
% figure, myMontagemn(asind(abs(mxypatptx2d)),5,2), caxis([0 1])
%figure, position_plots(asind(abs(mxypatptx2d)),[2 5],[0 fa],[],mask)

%save design4whole rf mxypatptx2d fa mask
%save design4partial rf grad mxypatptx2d fa mask targ
