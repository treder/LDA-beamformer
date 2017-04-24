# LDA-beamformer
Matlab code for LDA beamforming in EEG/MEG data

## Data-driven beamforming

Linear Discriminant Analysis (LDA) is a classifier that is popular for the decoding of event-related potentials (ERPs). In addition to this, [Treder et al (2016)](http://www.sciencedirect.com/science/article/pii/S1053811916000252) showed that LDA can also be used as a spatial filter to derive the time series of a brain source in EEG and MEG. In fact, LDA is mathematically equivalent to LCMV beamforming. Both methods differ only in the way that the spatial pattern is extracted (see figure).

![LDA vs beamforker](/images/LDA_vs_beamformer.jpg)

LDA beamformer and LCMV beamformer make use of a spatial pattern and the covariance matrix. In LCMV, the spatial pattern is derived from a source model, so it corresponds to the topography of a source voxel. In LDA beamforming, the spatial pattern is directly extracted from the EEG data (e.g. an ERP peak).

### Features of the LDA beamformer

* the LDA beamformer estimates the time series of a source
* it does not suffer from signal loss in the presence of correlated sources
* it can represent the time course of a distributed source (e.g. the common activity of both auditory cortices)
* it *cannot* be used for localising the source (you would use an ordinary beamformer for this)
* in simulations, [Treder et al (2016)](http://www.sciencedirect.com/science/article/pii/S1053811916000252) showed that it leads to better source estimates than the LCMV beamformer or PCA. The LCMV beamformer estimates time series in 2 steps (first localisation of the source voxel, then source estimation), whereas LDA beamformer skips the localisation step and estimates the time series directly.
* it can also be interpreted as a **one-class LDA** where the goal is to maximise the difference between a class and zero

## Using the LDA beamformer

The following gives a recipe about how to use the LDA beamformer. The code examples use Fieldtrip, but the examples can easily be adapted to other toolboxes such as EEGLAB. The LDA beamformer code itself is independent of any particular implementation.

**Step #1: extract spatial pattern from EEG/MEG data**

The LCMV beamformer obtains the spatial pattern of a source from a source model based on an MRI image. The LDA beamformer uses a spatial pattern that is derived from the EEG/MEG data itself. In the classical application of LDA as a classifier, the pattern represents a *difference* between two experimental conditions (pattern1 - pattern2). However, in LDA beamforming, a single experimental condition is sufficient. Here, it is simply the difference between the pattern and a *zero class* (pattern1 - 0) that is maximised. 

So how do we actually obtain the spatial pattern from the data? The following list gives some possible approaches:

* mean amplitude in an interval around an ERP peak (spatial pattern = ERP topography)
* a principal component of the covariance matrix
* a principal component of the (real part of the) cross-spectrum at a particular frequency

We will follow the first approach and obtain the spatial pattern from an interval around an ERP peak. Assuming the data is already epoched, the following code can be used for extracting the pattern in FieldTrip.

```Matlab
# Calculate ERP
erp = ft_timelockanalysis([], epoched);   

# We assume there is an ERP peak in the interval 0.3-0.5s
ival= [0.3 0.5];  

% To extract the spatial pattern and store it in a vector, we first find
% the corresponding time samples in the 0.3-0.5s window
time_idx = (erp.time >= ival(1) & erp.time <= ival(2));

% The mean across these time samples is our spatial pattern P
P = mean( erp.avg(:, time_idx),2);

```


**Step #2: obtain source time series**

We now want to obtain the time series of the (distributed) source associated with the pattern P that we just extracted. As inputs, the function only requires the spatial pattern P and continuous or epoched data. The data can come as a 2D (continuous) or 3D (epoched) matrix. Alternatively, a FieldTrip struct can be used. In this case, it is assumed that the EEG data is stored in the .trial field.

```Matlab
[w,dat_source,C] = LDAbeamformer(P,epoched);
```

The outputs are the spatial filter *w*, the source time series data *dat_source*, and the covariance matrix *C*.

**Step #3: analysis of source time series**

You can now use the source time series for any further analyses that are relevant to your experimental question. Examples:

* single-trial analysis (eg peak latencies/amplitudes)
* connectivity analysis (eg between the source time series and another source)
* multi-modal analysis (eg correlation between source time series and another concurrently recorded modality like fMRI or MEG)
* time-frequency analysis

## Nulling constraints

The beamformer tries to suppress all sources other than the source corresponding to the target spatial pattern. In some cases, we are interested in the relationship between two (or more) different sources, e.g. by way of their connectivity. In these cases, it is important to assure that the activity from source 1 does not leak into the estimate of source 2, and vice versa. This can be realised using nulling constraints (see eg [Nulling beamformer](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2818446/)). A nulling constraint is associated with additional sources. The activity from these sources is fully suppressed.

The LDA beamformer code also implements nulling constraints. If the input pattern P consists of a matrix instead of a vector, then the first pattern in this matrix (the first column) is considered as the target pattern and all the other patterns are "nulled" out, i.e., their activity is fully suppressed. In the following example, we are interested in two sources associated with two spatial patterns *P_source1* and *P_source2*. The second source is nulled out by providing it as an extra input pattern as follows:

```Matlab
[~,dat_source1] = LDAbeamformer([P_source1, P_source2],epoched); % activity of source 1 with source 2 fully suppressed
[~,dat_source2] = LDAbeamformer([P_source2, P_source1],epoched); % activity of source 2 with source 1 fully suppressed
```

