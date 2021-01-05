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

2. Process your experiment. This does three things: (1) loads the `block` and `Timeline` structs into `expInfo`, (2) generates `behavioralData` to hold trial-by-trial wheel movements and event times, and (3) generates `neuralData` to hold full traces as well as trial-by-trial event-aligned responses for each cell
```matlab
[expInfo, neuralData, behavioralData] = processExperiment(expInfo);
```

3. Retrieve data from video ROIs that were processed in Facemap
```matlab
eyeData = getEyeData(expInfo);
```
_NB: The identity/indexing of video ROIs will depend on how you decide to process your videos in Facemap. The conventions for toupee are described below_

This code uses cross-correlation to synchronize video frame times and global Timeline times, by comparing the SVD motion energy of an ROI containing the lick spout and optical beam events that were logged at the lick spout itself. An output plot lets you inspect the sync by overlaying video motion energy and optical lick events. 

We use this analysis since UDP signal times were not logged between eye & experiment computers at the B2 imaging rig. If you have UDP times available in your `Timeline` struct, these can be used to sync the two time series instead. 

4. Generate trial-by-trial event-aligned activity for each video ROI and store it in `eyeData`
```matlab
[eyeData] = alignFace(expInfo, eyeData, behavioralData);
```
# 2. Indexing trial types
Most data analyses will require comparing some types of trials to other types of trials (e.g., correct _vs._ incorrect, left _vs._ right). Generally, the repository calls these trial _conditions_ and there are a few scripts that will let you index the trials of your choosing.

## `contrasts = getUniqueContrasts(expInfo)` is a helper function that gives you a 1xn vector of all contrasts used in the session

## `initTrialConditions` lets you call specific name-value pairs to identify the exact conditions you want to focus on. 
Available pairs are:
  * `'repeatType'`: `{'all'}`, `{'random'}`, or `{'baited'}`
  * `'movementDir'`: `{'all'}`, `{'cw'}`, or `{'ccw'}` ('cw' refers to the movement a mouse would make to correctly report a left-side stimulus)
  * `'movementTime'`: `{'all'}`, `{'early'}`, or `{'late'}` (refers to when the mouse made its first movement with respect to the cue delay)
  * `'highRewardSide'`: `{'all'}`, `{'left'}`, or `{'right'}` (refers to which stimulus side had a high-value reward when reported correctly)
  * `'responseType'`: `{'all'}`, `{'correct'}`, or `{'incorrect'}`
  * `'rewardOutcome'`: `{'all'}`, `{'rewarded'}`, or `{'unrewarded'}` (this is useful for a 2AUFC task or a task where the reward valve fires probabilistically)
  * `'pastStimulus'`: `{'all'}`, `{'left'}`, `{'right'}`, or `{'zero'}` (refers to stimulus side)
  * `'pastMovementDir'`: `{'all'}`, `{'cw'}`, or `{'ccw'}`
  * `'pastResponseType'`: `{'all'}`, `{'correct'}`, or `{'incorrect'}`
  * `'trialsBack'`: integer (used in conjunction with 'past' conditions; default = 0)
  * `'switchBlocks'`: `{'all'}`, `{'beforeLeft'}`, `{'beforeRight'}`, `{'afterLeft'}`, or `{'afterRight'}` (selects the last (first) 50 trials before (after) a switch to a different reward block)
  * `'whichTrials'`: indexing vector or `{'all'}` (lets you select a custom range of trials)
  * `'specificRTs'`: 1 x 2 vector ([min max]) or `{'all'}` (selects trials that fall within the range you specify)

By default, running `trialConditions = initTrialConditions()` chooses 'all' for each trial condition. Name-value pairs can be concatenated with a semicolon and can be listed in any order.

Example: `trialConditions = initTrialConditions('responseType,{'correct'}; 'movementDir',{'cw'})` selects trials where the mouse was correct AND moved the wheel clockwise.

Selecting contrasts is done in a separate step.
  
3. 

#
# Conventions


