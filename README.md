# toupee
This toolbox is for processing and analyzing 2P imaging data aquired during behavior. It should give you everything you need to analyze how behavioral events are related to your imaging data.

# 0. What you'll need:
1. neural data: processed .mat files from [Suite2P](https://github.com/MouseLand/suite2p)
2. behavioral event data: a block.mat file from [Rigbox](https://github.com/cortex-lab/Rigbox)
3. experimental metadata: a Timeline.mat file from [Rigbox](https://github.com/cortex-lab/Rigbox)
4. video data: an eye_proc.mat file from [facemap](https://github.com/MouseLand/FaceMap) (if you want to analyze videos)

# 1. Load and process data
1. Initialize your experimental session by creating the `expInfo` struct. This will hold most of the experimental metadata (e.g., stimulus contrasts & onset times) and some behavioral data (e.g., trial-by-trial choices)
```matlab
expInfo = initExpInfo({{'LEW031'}},{{'2020-02-03',1,[1]}});
```
_BETA: there is functionality to processing multiple and/or cell-matched sessions, but this is still in development_
2. Process your experiment. This does three things: (1) loads the `block` and `Timeline` structs into `expInfo`, (2) generates `behavioralData` to hold trial-by-trial wheel movements and event times, and (3) generates `neuralData` to hold full traces as well as trial-by-trial event-aligned responses for each cell
```matlab
[expInfo, neuralData, behavioralData] = processExperiment(expInfo);
```
_BETA: there is functionality to processing multiple and/or cell-matched sessions, but this is still in development_
3. Retrieve data from video ROIs that were processed in Facemap
```matlab
eyeData = getEyeData(expInfo);
```
_NB1: The identity/indexing of video ROIs will depend on how you decide to process your videos in Facemap. The conventions for toupee are described below_

This code uses cross-correlation to synchronize video frame times and global Timeline times, by comparing the SVD motion energy of an ROI containing the lick spout and optical beam events that were logged at the lick spout itself. An output plot lets you inspect the sync by overlaying video motion energy and optical lick events. 

We use this analysis since UDP signal times were not logged between eye & experiment computers at the B2 imaging rig. If you have UDP times available in your `Timeline` struct, these can be used to sync the two time series instead. 

4. Generate trial-by-trial event-aligned activity for each video ROI and store it in `eyeData`
```matlab
[eyeData] = alignFace(expInfo, eyeData, behavioralData);
```
# 


