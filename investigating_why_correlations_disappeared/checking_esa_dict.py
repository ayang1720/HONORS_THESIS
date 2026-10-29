'''
does the esa pickle file with the 100 smois members actually have 100 smois members?
'''

import pickle
import numpy as np
esa_dict=pickle.load(open('ESA_with_100_members.pkl','rb'))
for i in esa_dict:
    print(i)
    print(np.shape(esa_dict[i]))