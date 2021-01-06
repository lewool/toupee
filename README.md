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
# 2. Plotting behavior
A few scripts for plotting the choice performance over a session are available.

#### `plotPsychometric(varargin)`
Plots a psychometric curve for a 2AFC task (one curve) or a biased-block 2AFC task (two curves). It reads the block file to automatically determine which task was run. 
A single session can be plotted by calling either:
  * `plotPsychometric({{'LEW025'}}, {{'2019-11-15',1}})`
  * `plotPsychometric(expInfo)`
  
Multiple sessions are concatenated by calling `plotPsychometric(mouseList, expList)`, where:
  * `mouseList = {{'Mouse1'}} or {{'Mouse1'},{'Mouse2'},{MouseN'}}`
  * `expList = {{'2018-06-10',2,[2 3]},{'2019-03-27',1,[1]}}` (NB: the final vector typically holds the same value as the second integer. There are occasional old files that concatnate over multiple sessions (hence the vector), but these are rarely encountered)

Multiple sessions can also be plotted by calling `plotPsychometric(expInfo)`, if `expInfo` is a struct of multiple experiments (see `initExpInfo.m`)

In the plot, green curves correspond to high-value left blocks, and orange curves correspond to high-value right blocks.

#### `pRight = plotBlockShifts(varargin)`
Plots a rolling mean of the proportion of right choices over time. This is useful when inspecting a biased-block 2AFC where you expect the mouse to go through epochs of more left or more right choices. It also outputs `pRight`, which is a vector of the rolling mean of p(right) over time.

In the plot, green corresponds to high-value left blocks, and orange corresponds to high-value right blocks.

This function can be called the same way as `plotPsychometric(varargin)`, above.

# 3. Indexing trial types
Most data analyses will require comparing some types of trials to other types of trials (e.g., correct _vs._ incorrect, left _vs._ right). Generally, the repository calls these trial _conditions_ and there are a few scripts that will let you index the trials of your choosing.

#### `trialConditions = initTrialConditions(names, values)` 
Prepares specific name-value pairs to identify the trial conditions you want to include in your analysis. Available pairs are:
  * `'repeatType'`: `'all'`, `'random'`, or `'baited'` ('baited' trials are those that were repeated after an incorrect response)
  * `'movementDir'`: `'all'`, `'cw'`, or `'ccw'` ('cw' refers to the movement a mouse would make to correctly report a left-side stimulus)
  * `'movementTime'`: `'all'`, `'early'`, or `'late'` (refers to when the mouse made its first movement with respect to the cue delay)
  * `'highRewardSide'`: `'all'`, `'left'`, or `'right'` (refers to which stimulus side had a high-value reward when reported correctly)
  * `'responseType'`: `'all'`, `'correct'`, or `'incorrect'`
  * `'rewardOutcome'`: `'all'`, `'rewarded'`, or `'unrewarded'` (this is useful for a 2AUFC task or a task where the reward valve fires probabilistically)
  * `'pastStimulus'`: `'all'`, `'left'`, `'right'`, or `'zero'` (refers to stimulus side)
  * `'pastMovementDir'`: `'all'`, `'cw'`, or `'ccw'`
  * `'pastResponseType'`: `'all'`, `'correct'`, or `'incorrect'`
  * `'trialsBack'`: integer (used in conjunction with 'past' conditions; default = 0)
  * `'switchBlocks'`: `'all'}`, `'beforeLeft'}`, `'beforeRight'`, `'afterLeft'`, or `'afterRight'` (selects the last (first) 50 trials before (after) a switch to a different reward block)
  * `'whichTrials'`: indexing vector or `'all'` (lets you select a custom range of trials)
  * `'specificRTs'`: 1 x 2 vector ([min max]) or `'all'` (selects trials that fall within the range you specify)

By default, running `trialConditions = initTrialConditions()` chooses 'all' for each trial condition. Name-value pairs can be concatenated with a semicolon and can be listed in any order.

Example: `trialConditions = initTrialConditions('responseType,'correct'; 'movementDir','cw')` selects trials where the mouse was correct AND moved the wheel clockwise.

Specifying contrast conditions is called outside of `trialConditions`; see below.
  
#### `[~, trialIDs] = selectCondition(expInfo, contrasts, behavioralData, trialConditions)` 
Generates a 1 x n vector of trial numbers that pass your conditions.

`contrasts` can be called as an integer or a vector, but must include a value that was used in the task. To check which contrasts were used in the task, run the helper function `contrasts = getUniqueContrasts(expInfo)`, which generates a vector of all contrasts used in the current session

Example: `[~, trialIDs] = selectCondition(expInfo, [-1 1], behavioralData, trialConditions)` generates the trial IDs for all trials with contrast = 1 or -1 AND passing any conditions you specified earlier in `trialConditions`

#### `trialTypes = getTrialTypes(expInfo, behavioralData, movementTime)`
A quick way to generate lots of different groups of conditions for use later on

Trials are automatically segregated by contrast, movementDir, responseType, and highRewardSide'. The only extra condition you can specify from the command line is `movementTime` ('early', 'late', or 'all'). Trial history is not supported.

Trials are organized by single conditions, interacting conditions, and contrast-corrected interacting conditions (the number of trials in each 'bin' is equalized, e.g., same number of correct vs incorrect -100% contrast trials)

# Conventions
Contrasts: Stimulus contrast is expressed on a scale from 0–1, and is signed to denote the screen it appeared on. Left = negative; right = positive.

Choices: Animal choices are either -1 ('chose left') or +1 ('chose right').'Chose left' means that the mouse reported a stimulus on the left by turning the wheel CW. 'Chose right' means the mouse reported a stimulus on the right by turning the wheel CCW. 

CW vs CCW: These are wheel directions, determined from the perspective of the mouse. CW is the wheel action that moves a stimulus to the right; CCW is the wheen action that moves a stimulus to the left. These designations can be used independently of correct/incorrect. CW turns _increase_ the rotary encoder value; CCW turns _decrease_ the rotary encoder value (raw encoder values aren't used much but it's good to bear in mind).

Colors:
  * Green vs orange: This color pair is used to compare blocks of high-value left choices (green) to blocks of high-value right choices (orange). 
  * Blue vs red: This color pair is used to compare stimulus position or brain hemisphere. Blue  means left (or contralateral); red means right (or ipsilateral).
  * Green vs brown: This color pair is used to denote correct (green) versus incorrect (brown) trial outcomes.

