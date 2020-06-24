# +toupee

Contains the toupee repository's source code for organizing, analyzing, and plotting behavioral and neural data.

The base datatype of this repository is the `expInfo` struct array returned by `+meta\processExperiment.m`. This struct array contains a struct for each experiment session, where each struct contains meta, behavioral, and neural data for that session. The behavioral and neural data are respectively stored in the `behavioralData` and `neuralData` sub structs. The possible fields in the `behavioralData` struct are `trials` (filtered trials), `eventTimes` (signals and daq times for events of interest),... The possible fields in the `neuralData` struct are... 

Note on function signatures: All the `+meta\` and 'get' functions in `+behavioral\` and `+neural\` require an `expInfo` struct input (though they can output arrays in addition to returning an updated `expInfo` struct). The `+plot\` and analysis functions  in `+behavioral\` and `+neural` accept and return base matlab datatypes (typically cell arrays and/or arrays).

Contents:
`+behavioral\` : Contains code for organizing and computing stats on behavioral data.
`+meta\` : Contains code for returning meta-information about experiments.
`+misc\` : Contains miscellaneous helper functions.
`+neural\` : Contains code for organizing and computing stats on neural data ([suite2p output](https://github.com/MouseLand/suite2p)).
`+plot\` : Contains code for plotting behavioral and neural data.