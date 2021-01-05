# toupee
This toolbox is for processing and analyzing 2P imaging data aquired during behavior. It should give you everything you need to analyze how behavioral events are related to your imaging data.

# 0. What you'll need:
1. neural data: processed .mat files from [Suite2P](https://github.com/MouseLand/suite2p)
2. behavioral event data: a block.mat file from [Rigbox](https://github.com/cortex-lab/Rigbox)
3. experimental metadata: a Timeline.mat file from [Rigbox](https://github.com/cortex-lab/Rigbox)
4. video data: an eye_proc.mat file from [facemap](https://github.com/MouseLand/FaceMap) (if you want to analyze videos)

# 1. Loading and preprocessing data
1. Initialize your experimental session by creating the `expInfo` struct. This will hold most of the experimental metadata (e.g., stimulus contrasts & onset times) and some behavioral data (e.g., trial-by-trial choices)
```matlab
expInfo = initExpInfo({{'LEW031'}},{{'2020-02-03',1,[1]}});
```
# 
