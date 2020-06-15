# toupee

This toolbox is for processing and analyzing neural data (fluorescence traces from calcium imaging) and behavioral data (in which the subject manipulates a rotary encoder to indicate response) acquired from alternative forced choice tasks. It should give you everything you need to analyze how behavioral events are related to your imaging data.

# What you'll need:

1. processed .mat files from [Suite2P](https://github.com/cortex-lab/Suite2P)
2. a block.mat file from [Rigbox](https://github.com/cortex-lab/Rigbox)
3. a Timeline.mat file from [Rigbox](https://github.com/cortex-lab/Rigbox)

# What's in here

## Repository structure:

```
+toupee\
      +behavior\
      +data\
      +meta\
      +misc\
      +neural\
archive\
docs\
examples\
tests\
```

The source code that processes and analyzes data is organized in [packages](https://uk.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html) within `+toupee\`. `archive\` contains old source code and examples, `docs\` contains files that detail meta-information about this repository, `examples\` contains files that show usage examples of the source code, and `tests\` contains unit and integration tests for the source code.

Each directory in this repository contains its own README file summarizing that directory's purpose and contents.

# Getting Started

First, rename the `+toupee\+data\getPathsTemplate.m` function to `+toupee\+data\getPaths.m`, and then change the code within this function to return paths that point to your subject data: each path returned by this function should contain the following nested directory structure: `\subject\data\session` (for example, the first path in `getPathsTemplate.m` is `\\znas.cortexlab.net\Subjects`, which contains data for LEW031's first session on 2020-02-03 in `\\znas.cortexlab.net\Subjects\LEW031\2020-02-03\1`). 

After setting your paths, open and follow the instructions in the `behavior_analysis`and `neural_analysis` scripts in `examples\` to learn how to get started using the code in this repository.