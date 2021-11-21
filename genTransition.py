import sys
import numpy as np
import pandas as pd
df=pd.read_csv(sys.argv[1]+'.csv')
print(df);
P_mat0=np.array([[0]*8]*8);
P_mat=np.array([[0]*8]*8,dtype=np.float16);
for i in range(0,len(df)-1):
    P_mat0[int(df.iloc[i,2])-1,int(df.iloc[i+1,2])-1]+=1;

for j in range(0,7):
    for i in range(0,7):
        P_mat[i,j]=P_mat0[i,j]/(P_mat0.sum(axis=1)[j]);
print(P_mat);
pd.DataFrame(P_mat).to_csv("./"+sys.argv[1]+"_transitionMat.csv")
