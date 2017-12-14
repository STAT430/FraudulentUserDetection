

#python victim.py test-1k-1k-random.txt post-1k-1k-random-more-iterations.txt &
#python victim.py testfile.txt scorefile.txt trainfile.txt


if __name__ == "__main__":
    import sys
    import math
    import os
    import operator
    import numpy as np

    
    test_file = open(sys.argv[1])
    actives = set([])
    suspends = set([])
    
    line = test_file.readline()
    line_parts = line[:-1].split(' ')
    for id in line_parts:
        try:
            if len(id) != 0:
                actives.add(int(float(id)))
        except Exception:
            pass

    line = test_file.readline()
    line_parts = line[:-1].split(' ')
    for id in line_parts:
        if len(id) != 0:
            suspends.add(int(float(id)))
#    print len(actives)
#    print len(suspends)

    scores = {} # Dict
    score_file = open(sys.argv[2])
    
    for line in score_file:
        line_parts = line.split(' ')
        node = int(float(line_parts[0]))
        s = float(line_parts[1])
        #s = 1 - float(line_parts[1])

        if node in actives or node in suspends:
            scores[node] = s
    
    # Default ascendent 
    sorted_scores = sorted(scores.items(), key=operator.itemgetter(1)) 
    
    
    #for k in [1000, 10000, 50000, 100000, 1000000, 10000000]:
#    for k in [10, 30, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]:
    #for k in [1, 2, 3, 4] :
    #interval = 20000
    interval = int(sys.argv[3])
    num = int(sys.argv[4])

    n = 0
    for k in np.array(range(1,num+1))*interval:
    #for k in [interval, 2*interval, 3*interval, 4*interval, 5*interval, 6*interval, 7*interval, 8*interval, 9*interval, 10*interval]:
        prev_n = n
        if k>len(sorted_scores):
            continue
        
        index_low = k - 1
        index_high = k - 1 
        
        i = k - 1
        while i >= 1 and sorted_scores[i - 1][1] == sorted_scores[i][1]:
            i -= 1
        index_low = i
        
        
        i = k - 1
        while i < len(sorted_scores) - 1 and sorted_scores[i][1] == sorted_scores[i + 1][1]:
            i += 1
        index_high = i 
        
        num_sybil = 0.0
        for pair in sorted_scores[index_low: index_high + 1]:
            if pair[0] in suspends:
                num_sybil += 1.0
        
        n = num_sybil / (index_high - index_low + 1.0) * (k - index_low)
        
        for pair in sorted_scores[:index_low]:
            if pair[0] in suspends:
                n += 1
    
        
        print(interval, (n-prev_n))

#    print ''

#    #for k in [1000, 10000, 50000, 100000, 1000000, 10000000]:
#    for k in [10, 30, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]:
#    #for k in [400]:
#        if k>len(sorted_scores):
#            continue
#
#        index_low = -k
#        index_high = -k
#           
#        i = -k
#        while i >= -len(sorted_scores) + 1 and sorted_scores[i - 1][1] == sorted_scores[i][1]:
#            i -= 1
#        index_low = i
#           
#           
#        i = -k
#        while i < - 1 and sorted_scores[i][1] == sorted_scores[i + 1][1]:
#            i += 1
#        index_high = i 
#           
#        num_sybil = 0.0
#        for pair in sorted_scores[index_low: index_high]:
#            if pair[0] in suspends:
#                num_sybil += 1.0
#        if sorted_scores[index_high][0] in suspends:
#                num_sybil += 1.0
#                   
#        n = num_sybil / (index_high - index_low + 1.0) * (index_high + k + 1)
#           
#        if index_high != -1:
#            for pair in sorted_scores[index_high + 1:]:
#                if pair[0] in suspends:
#                    n += 1
#           
#        print k, n

        
        
    
    
    
