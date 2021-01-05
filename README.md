# toupee
This toolbox is for processing and analyzing 2P imaging data aquired during behavior. It should give you everything you need to analyze how behavioral events are related to your imaging data.

# What you'll need:
1. processed .mat files from [Suite2P](https://github.com/MouseLand/suite2p)
2. a block.mat file from [Rigbox](https://github.com/cortex-lab/Rigbox)
3. a Timeline.mat file from [Rigbox](https://github.com/cortex-lab/Rigbox)
4. an eye_proc.mat file from [facemap](https://github.com/MouseLand/FaceMap) (if you want to analyze videos)

# What's in here
- /behavior
  - loadData.m: loads the block.m and 
  - selectCondition.m: chooses the relevant trials based on a list of qualifiers
  - getEventTimes.m: extracts behavioral events from Timeline in order to compare to imaging events
  - wheelAnalysis: a toolbox for extracting wheel movements (copied from [here](https://github.com/cortex-lab/wheelAnalysis) and probably out of date)
- loadExpTraces.m: load the Suite2P data structures
- getExpTraces.m: chop up & align the traces to a predefined behavioral event
- getPlaneFrameTimes.m: assign a Timeline time to each frame acquired during your imaging experiment
- chooseCellType.m: pick cells to look at based on rough-&-tumble criteria (visual, movement, or all)
- plotResps_B2AFC.m: plots cell responses, divided by visual contrast & reward block (grand mean & browser)
# 
