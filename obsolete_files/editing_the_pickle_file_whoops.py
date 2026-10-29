'''
6/18/26
accidentally made the dictionary in avg_dict size 40962 when it should've been 729
this py file corrects the mistake by removing the 40233 superfluous elements in avg_dict's array
'''

import pickle

#grabbing data
avg_dict=pickle.load(open('/fs/scratch/PAS3252/yang/HONORS_THESIS/averages_100.pkl','rb'))
for i in avg_dict:
    print(len(avg_dict[i]))
    avg_dict[i]=avg_dict[i][:729]
pickle.dump(avg_dict,open('/fs/scratch/PAS3252/yang/HONORS_THESIS/averages_100.pkl','wb'))