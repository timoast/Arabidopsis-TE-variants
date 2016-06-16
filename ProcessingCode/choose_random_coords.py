import random
import sys


chroms = {'chr1': 30427671, 'chr2': 19698289,'chr3':23459830,'chr4':18585056,'chr5':26975502}

def genRandCoords(chroms, number, l):
    for x in range(int(number)):
        c = random.choice(chroms.keys())
        start = random.randrange(chroms[c]-int(l))
        end = start + int(l)
        print(c+'\t'+str(start)+'\t'+str(end))

genRandCoords(chroms, sys.argv[1], sys.argv[2])