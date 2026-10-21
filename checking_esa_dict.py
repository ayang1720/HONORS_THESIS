import pickle
import numpy as np
esa_dict=pickle.load(open('ESA_with_100_members.pkl','rb'))
for i in esa_dict:
    print(i)
    print(np.shape(esa_dict[i]))