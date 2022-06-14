# LOHC-Score
A dehydrogenation scoring system for LOHC candidates based on 2D molecular features


## Installation
Clone or download score.py (found in this repo), and then install it's dependencies:

This script uses three [AutoMech](https://tcg.cse.anl.gov/papr/codes/automech.html) libraries: [autochem](https://sne-autochem.readthedocs.io/en/latest/), [autoio](https://sne-autoio.readthedocs.io/en/latest/) and [autofile](https://sne-autofile.readthedocs.io/en/dev/).  The libraries can be installed by cloning/forking the [AutoMech Github](https://github.com/Auto-Mech) or through the [AutoMech conda](https://anaconda.org/auto-mech): 
```bash
>>> conda env create auto-mech/amech-env
>>> conda activate amech-env
```

## Usage
The scoring system script requires an input csv file of InChI strings with tanimoto similarity scores (or any random float) to be prepared beforehand.  Then it is a simple command line execution.
```
>>> conda activate amech-env
>>> python score.py -i input.csv -o output.csv
```

```
>>> cat input.csv
"InChI=1S/C13H9NS2/c15-13-12(10-6-2-1-3-7-10)16-11-8-4-5-9-14(11)13/h1-9H/p+1",5.3,0.35
"InChI=1S/C12H6N2/c1-2-4-10-9(3-1)13-11-6-7-5-8(7)12(11)14-10/h1-6H",4.5,0.366667
"InChI=1S/C12H11N/c1-2-3-9-13-10-8-11-6-4-5-7-12(11)13/h1,4-8,10H,3,9H2",4.5,0.4
"InChI=1S/C12H11N/c1-2-3-7-11-9-10-6-4-5-8-12(10)13-11/h1,4-6,8-9,13H,3,7H2",4.5,0.428571
"InChI=1S/C12H9N/c1-2-3-7-11-9-10-6-4-5-8-12(10)13-11/h2,4-6,8-9,13H,1H2",4.5,0.428571

>>> cat output.csv
InChI				, factor,    score,tanimoto,% wt h2,sp3_carb,sp2_carb,sp_carb,M-no5mem,  M-5mem,M-nitcont, M-SOcont,B-no5mem,  B-5mem,B-nitcont, B-SOcont,P-no5mem,  P-5mem,P-nitcont, P-SOcont,nonCfuse,SObadpos,SO_adjrng,3mem_bi+,one-posi,onethree
"InChI=1S/C12H6N2/c1-2-4-10-9(3-1)13-11-6-7-5-8(7)12(11)14-10/h1-6H",9.532,4.5,7.01609,0,12,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,4,6
"InChI=1S/C12H9N/c1-2-3-7-11-9-10-6-4-5-8-12(10)13-11/h2,4-6,8-9,13H,1H2",9.157,4.5,7.02703,0,10,2,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,2,2
"InChI=1S/C12H11N/c1-2-3-9-13-10-8-11-6-4-5-7-12(11)13/h1,4-8,10H,3,9H2",9.075,4.5,5.92525,2,8,2,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,2,1
"InChI=1S/C12H11N/c1-2-3-7-11-9-10-6-4-5-8-12(10)13-11/h1,4-6,8-9,13H,3,7H2",8.813,4.5,5.92525,2,8,2,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,2,2
```

## Features
Our current descriptors and their additive value, as of June 9

- Cutoff threshold of 5.5 %wtH2 is treated as a cutoff parameter rather than a descriptor
- % wt h2: 0, this feature is treated as a cutoff parameter so any under 5.5 get total scores of 0
- sp3_carb: 0.1153, number of sp3 carbons
- sp2_carb: 0.2872, number of sp2 carbons
- sp_carb:  1.1298, number of sp carbons
- M-no5mem: -1.2127, number of monocyclic rings that aren't 5-membered rings
- M-5mem:   -1.6513, number of 5-membered monocyclic rings
- M-nitcont: -0.0955, number of monocyclic rings with a nitrogen atom in the ring
- M-SOcont: 0.9387, number of monocyclic rings with a sulfur or oxygen atoms in the ring
- B-no5mem: 5.2805, number of bicyclic rings that don't have any 5-membered rings
- B-5mem:   5.0358, number of 5-membered ring containing bicyclic rings
- B-nitcont: -.3393, number of biyclic rings with a nitrogen atom in any ring
- B-SOcont: -.0191, number of bicyclic rings with a sulfur or oxygen atoms in any ring
- P-no5mem: 1.0000, number of polycyclic (3 or more) rings that don't have any 5-membered rings
- P-5mem:   3.4910, number of 5-membered ring containing polycyclic rings
- P-nitcont: 0.5179, number of polyyclic rings with a nitrogen atom in any ring
- P-SOcont: .06033, number of polycyclic rings with a sulfur or oxygen atoms in any ring
- nonCfuse, 0.5200, number of non-carbon atoms that are the fusion point of two rings
- SObadpos, 0.9095, number of times a sulfur or oxygen is gamma to a fusion point
- SO_adjrng, 2.000, number of times two fused rings both have a sulfur or oxygen
- 3mem_bi+, 3.9448, number of times a 3-membered ring is fused to another ring
- one-posi, -0.0735, number of substituents on a ring (including another ring)
- onethree, -0.2623, number of 1,3 substituent or 1,3 nitrogen interactions on a ring (including another ring) onethree]}
