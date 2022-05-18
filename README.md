## OED parameters
### Constant Flip angle

  1. Pyruvate FA = 35, Lactate FA = 28 [^1]  
  2. Pyruvate FA = 3 , Lactate FA = 28 [^2]  

### Variable Flip angle
  
  3. worktmpLB [^1] = load('poptNG5Nu3interior-pointTotalSignalSNR02Hermite.mat');
```
     >> worktmpLB.params.FaList
```
  4. worktmpUB [^2] = load('poptNG5Nu3interior-pointTotalSignalSNR20Hermite.mat');
```
     >> worktmpUB.params.FaList
```

### Model parameters 
```
     >> worktmpLB.params
```

[^1]: $\kk_{OED_{20}}$ 
[^2]: $\kk_{OED_{2}}$ 
