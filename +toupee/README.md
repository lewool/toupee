# +toupee

Contains the toupee repository's source code for organizing, analyzing, and plotting behavioral and neural data.

## Data Architecture

The base datatype of this repository is `expInfo`, a table returned by `+meta\processExperiment.m`. This table contains a row for each experiment session, and columns for meta, behavioral, and neural data for that session. 

Raw data saved into files from the experiment session (found in the session's directory) can be added to and accessed directly within the `expInfo` table: this data can be loaded via `+meta\loadDatafile.m`, which will make it available in `expInfo` in columns that end with the name `File`. Any number of raw data files can be loaded/added in such a manner.

Specific behavioral and neural data are respectively stored in `expInfo`'s `behavioralData` and `neuralData` sub tables. 

*Note on function signatures*: `+meta\loadDatafile.m` and the 'get' functions in `+behavioral\` and `+neural\` require an `expInfo` table input, and they can return base matlab datatypes in addition to returning an updated `expInfo` table. The `+plot\` and analysis functions  in `+behavioral\` and `+neural` accept and return base matlab datatypes.

## Contents:

`+behavioral\` : Contains code for organizing and computing stats on behavioral data.
`+meta\` : Contains code for returning meta-information about experiments.
`+misc\` : Contains miscellaneous helper functions.
`+neural\` : Contains code for organizing and computing stats on neural data ([suite2p output](https://github.com/MouseLand/suite2p)).
`+plot\` : Contains code for plotting behavioral and neural data.