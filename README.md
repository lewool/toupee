# toupee

This MATLAB toolbox is for processing and analyzing neural data (suite2p-processed fluorescence traces from calcium imaging) and behavioral data (Rigbox-processed sensor and actuator signals from a task in which a subject manipulates a rotary encoder to indicate response) acquired in behavioral neuroscience experiments. It should give you everything you need to analyze how behavioral events are related to your imaging data.

## Required data

1. neural data: processed .mat files from [Suite2P MATLAB](https://github.com/cortex-lab/Suite2P), or .npy files from [Suite2P Python](https://github.com/MouseLand/suite2p) for a given experiment session.

2. behavioral data: a block.mat and a Timeline.mat file from [Rigbox](https://github.com/cortex-lab/Rigbox) for a given experiment session.

*Note*: This repository assumes that all your subject data is contained in directories which hold the following nested directory structure: `\subject\date\session`. Using this directory structure, the code can find the experiment data for any session of any date of any subject across directories (e.g. you may have all the data for one subject saved on one server, and the data for another subject saved on another server, but as long as these servers contain the same nested directory structure, this code can search both).

## Required software

1. [MATLAB 2019b or later](https://www.mathworks.com/downloads/)
2. [Git](https://git-scm.com/downloads)

## Contents

### Repository structure:

```
+toupee\
      +behavioral\
      	+wheel\
      +meta\
      	+npy\
      +misc\
      +neural\
      +plot\
archive\
docs\
examples\
tests\
wheelAnalysis\
```

The source code that processes and analyzes data is organized in [packages](https://uk.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html) within `+toupee\`. As for the other folders, `archive\` contains old source code and examples, `docs\` contains files that detail meta-information about this repository, `examples\` contains files that show usage examples of the source code, and `tests\` contains unit and integration tests for the source code.

Each subdirectory in this repository contains its own `readme` file summarizing that directory's purpose and contents.

## Getting started

1. Clone the repository: Open your system terminal, navigate into the target directory in which you wish to clone this repository, and run the git command `git clone --recurse-submodules https://github.com/lewool/toupee`.

2. Set paths to your data: You must create a `getPaths.m` function within `+toupee\+data\` that contains the paths to your data: when searching for paths to your data, the code in the repository tries to find and call this function. You can use `+toupee\+data\getPathsTemplate.m` as a guide.

*Note*: Each path returned in your personal `getPaths.m` function should contain the following nested directory structure: `\subject\date\session` (for example, the first path in `getPathsTemplate.m` is `\\znas.cortexlab.net\Subjects`, which contains data from the first session on 2020-02-03 for subject LEW031 in `\\znas.cortexlab.net\Subjects\LEW031\2020-02-03\1`)

3. Run examples: In MATLAB, add `toupee\` and all its subfolders to your path. Then open and follow the instructions in the `behavioralAnalysis` and `neuralAnalysis` scripts in `examples\`.

## Notes on naming conventions:

- 'experiment session', 'experiment', and 'session' all refer to the same thing: a single, continuous block of time during which a subject was performing a behavioral task and behavioral data was recorded via Rigbox and (possibly) neural data was recorded via ScanImage and later processed via suite2p.

- 'expRef': An 'experiment reference' string, used by Rigbox to refer to a specific experiment session. It takes the form of '<expDate>_<expNum>_<subject>', e.g. '2020-02-03_1_LEW031'.

- A right choice/action indicates choosing the stimulus that appears on the right side: for a task in which the stimulus has to be moved to the center of the screen, this means turning the wheel to the left (i.e. a left wheel movement).