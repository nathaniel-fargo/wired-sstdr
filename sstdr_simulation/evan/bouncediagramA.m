G1=-1/3; G2=-1/3; G3 = -1/3; G6=-1/3;
G4=1;G5=1;
V1p=1
V1m=V1p*G2
V2p = V1m*G3
V2m=V2p*G2
V3p=V2m*G3
V4p=V1p*(1+G2)
V4m=V4p*G4
V5p=V4m*G6
V6m=V4m*(1+G6)
V6p=V6m*G3
V7p=V4p
V7m=V7p*G5
V8m=V7m*(1+G6)
V9p=V8m*G3
P0= 1
P200=V1m+V2p
P300 = V8m+V9p
P400=V2m+V6m+V3p*V6p
