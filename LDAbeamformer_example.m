% The LDA beamformer is used to obtain the TIME SERIES OF A SOURCE. To this
% end, the spatial pattern of the source needs to be estimated from the
% EEG/MEG data. 
%
% This example makes use of the FieldTrip toolbox for preprocessing and
% plotting the data.
% We will use an example MEG dataset from FieldTrip, which is
% used in the tutorial on event-related fields: 
% http://www.fieldtriptoolbox.org/tutorial/eventrelatedaveraging
% Make sure that the MEG file is in your working directory.
clear all
close all

%% (preprocessing taken from the FieldTrip ERF tutorial - FIC condition)
% find the interesting segments of data
cfg = [];                                           % empty configuration
cfg.dataset                 = 'Subject01.ds';       % name of CTF dataset  
cfg.trialdef.eventtype      = 'backpanel trigger';
cfg.trialdef.prestim        = 1;
cfg.trialdef.poststim       = 2;
cfg.trialdef.eventvalue     = 3;                    % trigger value for fully incongruent (FIC)
cfg = ft_definetrial(cfg);            

% remove the trials that have artifacts from the trl
cfg.trl([15, 36, 39, 42, 43, 49, 50, 81, 82, 84],:) = []; 

% preprocess the data
cfg.channel    = {'MEG'};%, '-MLP31', '-MLO12'};        % read all MEG channels except MLP31 and MLO12
cfg.demean     = 'yes';
cfg.baselinewindow  = [-0.2 0];
cfg.lpfilter   = 'yes';                              % apply lowpass filter
cfg.lpfreq     = 35;                                 % lowpass at 35 Hz.

epoched_FIC = ft_preprocessing(cfg);   

%%% Timelock analysis
cfg = [];
avg_FIC = ft_timelockanalysis(cfg, epoched_FIC);

%% Plot one MEG channel
figure
cfg.xlim        = [-0.2 1.0];
cfg.ylim        = [-1e-13 3e-13];
cfg.channel     = 'MLC24';
ft_singleplotER(cfg,avg_FIC);

%% Obtaining the spatial pattern of the ERP peak
% In the channel plot, we see that there is a peak in the interval
% 0.45-0.75s. Let's mark this interval.
ival= [0.45 , 0.75];
h= rectangle('Position',[ival(1),cfg.ylim(1),diff(ival),diff(cfg.ylim)],'FaceColor',[1 1 1]*0.7);
uistack(h,'bottom');

% Plot the SPATIAL PATTERN (= topography) for this interval. This is the
% spatial pattern that we will use in the LDA beamformer.
figure
cfg = [];
cfg.xlim = ival;
cfg.colorbar = 'yes';
ft_topoplotER(cfg,avg_FIC);

% To extract the spatial pattern and store it in a vector, we first find
% the corresponding time samples in the 0.45-0.75s window, and then 
% calculate the average spatial pattern across these time samples
time_idx = (avg_FIC.time >= ival(1) & avg_FIC.time <= ival(2));
P = mean( avg_FIC.avg(:, time_idx),2);

%% LDA beamforming ------------------------------------------------
% To obtain the LDA beamformer, we need the spatial pattern P and the
% corresponding continuous or epoched data (NOT the averaged data).
[w,lda,C] = LDAbeamformer(P,epoched_FIC);

%% Plot covariance matrix, spatial pattern P, and spatial filter w
% Get channel layout for plotting
cfg             = [];
cfg.layout      = 'CTF151.lay';
cfg.skipscale   = 'yes';
cfg.skipcomnt   = 'yes';
lay= ft_prepare_layout(cfg);

% Use the low-level plotting function for plotting arbitrary vectors
cfg_topo= {'mask' lay.mask, 'datmask' [], 'interplim' 'mask'};
cfg_lay = { 'point' 'no' 'box' 'no' 'label' 'no' };

figure
subplot(1,3,1), 
    imagesc(C),title('Covariance matrix')
subplot(1,3,2), 
    ft_plot_topo(lay.pos(:,1), lay.pos(:,2), P, cfg_topo{:});
    ft_plot_lay(lay, cfg_lay{:})
    title('Spatial pattern P')
subplot(1,3,3), 
    ft_plot_topo(lay.pos(:,1), lay.pos(:,2), w, cfg_topo{:});
    ft_plot_lay(lay, cfg_lay{:})
    title('Spatial filter w')

%% Source-level ERP
figure
cfg = [];
avg_LDA = ft_timelockanalysis(cfg, lda);

cfg= [];
cfg.xlim        = [-0.2 1.0];
cfg.channel     = 'MLC24';
ft_singleplotER(cfg,avg_FIC);

% Scale LDA source output to fit the channel plot
yl= ylim;
hold all
plot(avg_LDA.time, avg_LDA.avg * yl(2))

set(gca,'ylim',[yl(1) max(yl(2), max(avg_LDA.avg * yl(2)) )])

legend({cfg.channel, avg_LDA.label{1}})
title('Event-related field')

%% Investigating the effect of the regularisation parameter gamma on the spatial filter w
% gamma is the regularisation parameter set between 0= no regularisation
% and 1=full regularisation
figure
gammas= [0.0001:0.15:1];
for ii=1:numel(gammas)
    w = LDAbeamformer(P,epoched_FIC, gammas(ii));
    
    subplot(1,numel(gammas),ii)
    ft_plot_topo(lay.pos(:,1), lay.pos(:,2), w, cfg_topo{:});
    ft_plot_lay(lay, cfg_lay{:})
    title(sprintf('w (gamma=%1.3f)',gammas(ii)))
end