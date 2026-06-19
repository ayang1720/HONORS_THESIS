'''
6/18/26s
'''

import pickle

#grabbing data
avg_dict=pickle.load(open('/fs/scratch/PAS3252/yang/HONORS_THESIS/averages_100.pkl','rb'))
for i in avg_dict:
    print(len(avg_dict[i]))
    avg_dict[i]=avg_dict[i][:729]
pickle.dump(avg_dict,open('/fs/scratch/PAS3252/yang/HONORS_THESIS/averages_100.pkl','wb'))