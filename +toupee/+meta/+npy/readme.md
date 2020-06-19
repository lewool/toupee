This directory contains code to read and write NumPy's npy format (`.npy` files) in MATLAB. Based on this repository: https://github.com/kwikteam/npy-matlab

## Contents

`constructNPYheader` :
`readNPY` :
`readNPYheader` :
`writeNPY` :

## Usage example

```matlab
>> import toupee.meta.npy.*
>> a = rand(5,4,3);
>> writeNPY(a, 'a.npy');
>> b = readNPY('a.npy');
>> sum(a(:)) == sum(b(:))
ans =

    1
```