#!/usr/bin/python
'''
replace specific atoms with other atoms. Now It can only change one atom once time.
Last updata at 1/31/2016, by Jincheng Liu
'''

import sys,math,linecache,time
import numpy as np



def atom_colored():

    output = ''.join([fileName.split('.')[0], '-colored', '.', fileName.split('.')[1]])
    output_write = []
    for i in range(0, loop_N):
        for j in range(0,int(numOfAtomToBeChanged)+1):
            output_write.append(lines[i * (int(atomNum) + 2) + j])
            #open(output,'a').write(lines[i * (int(atomNum) + 1) + j])
        tmpl = lines[i * (int(atomNum) + 2) + int(numOfAtomToBeChanged) + 1].split(' ')
        while '' in tmpl:
            tmpl.remove('')
        #print(tmpl)
        for k in tmpl:
            if k == 'Au':
                output_write.append(targitAtom + ' ')
                #open(output, 'a').write('Pt' + '\t')
            else:
                output_write[-1] = output_write[-1] + (k + ' ')
        for j in range(int(numOfAtomToBeChanged) + 2, int(atomNum) + 2):
            output_write.append(lines[i * (int(atomNum) + 2) + j])
    list2str = ''.join(output_write)
    open(output, 'w').write(list2str)
    print("new coordinates have been wrote to",output)
    return None


def atom_to_center_time():
    output = ''.join([fileName.split('.')[0], '-colored', '.', fileName.split('.')[1]])
    with open(output,'r') as f:
        lines = f.readlines()
        print('total lines of colored file:',len(lines))
    distance_mc_pt = []
    for i in range(0, loop_N):
        coord_sum = [0, 0, 0]
        count_Au = 0
        for j in range(0, int(atomNum) + 2):
            tmpl = lines[i * (int(atomNum) + 2) + j].split(' ')
            while '' in tmpl:
                tmpl.remove('')
            if tmpl[0] == targitAtom:
                tmpl.remove(targitAtom)
                coord_pt = list(map(lambda x: float(x), tmpl))
            if tmpl[0] == 'Au':
                count_Au += 1
                coord_sum[0] += float(tmpl[1])
                coord_sum[1] += float(tmpl[2])
                coord_sum[2] += float(tmpl[3])

        coord_mass_center = list(map(lambda x: x/count_Au, coord_sum))
        mc = np.array(tuple(coord_mass_center))
        pt = np.array(tuple(coord_pt))
        distance_mc_pt.append(str(np.linalg.norm(mc - pt)))
    print(len(distance_mc_pt))
    list2str = '\n'.join(distance_mc_pt)
    open('distance_mc_pt.txt', 'w').write(list2str)
    print('The distance between specific atom and cluster\'s mass center\
     are written to distance_mc_pt.txt')

def center_metal_distribution(interval = 0.1):
    with open(fileName,'r') as f:
        lines = f.readlines()
        #print('total lines:',len(lines))

    distance_mc_pt = []

    for i in range(0, int(loop_N/2.)):
        coord_sum = [0, 0, 0]
        count_Me = 0
        coord_Me = []
        for j in range(0, int(atomNum) + 2):
            tmpl = lines[i * (int(atomNum) + 2) + j].split(' ')
            while '' in tmpl:
                tmpl.remove('')
            if  tmpl[0] == 'Au':
                tmpl.remove('Au')
                count_Me += 1
                coord_sum[0] += float(tmpl[0])
                coord_sum[1] += float(tmpl[1])
                coord_sum[2] += float(tmpl[2])

                coord_Me.append(list(map(lambda x: float(x), tmpl)))
        coord_mass_center = list(map(lambda x: x/count_Me, coord_sum))

        for k in coord_Me:
            mc = np.array(tuple(coord_mass_center))
            me = np.array(tuple(k))
            distance_mc_pt.append(str(np.linalg.norm(mc - me)))

    distribution = list(map( lambda x: x*0 , [i for i in range(0, 150)]))
    for i in range(len(distribution)):
        for j in distance_mc_pt:
            if i * interval <= float(j) < (i+1) * interval:
                distribution[i] += 1

    # do not do number counting
    #list2str = '\n'.join(distance_mc_pt)
    #open('distribution_mc_pt.txt', 'w').write(list2str)

    list2str = '\n'.join(list(map(lambda x:str(x/sum(distribution)), distribution)))
    open('distribution_mc_pt.txt', 'w').write(list2str)
    print('The probability distribution functions P(rcm) \
    of the all Au atoms relative to the centre-of-mass \
    are written to distribution_mc_pt.txt')
    return None





print('#######  START  #######')
do_atom_colored = True
do_atom_to_center_time = True
do_center_metal_distribution = True
fileName = 'newpos.xyz'       # input filename
# filename = sys.argv[1]
numOfAtomToBeChanged = '233'  # order Of Atom To Be Changed
# numOfAtom = sys.argv[2]
targitAtom = 'Pt'       # the targit atom to become
# targitAtom = sys.argv[3]



atomNum = (linecache.getline(fileName,1)).strip()
print('total atom number:', atomNum)
print('the atom order to be changed:', numOfAtomToBeChanged, targitAtom)

#count the iteration number
with open(fileName,'r') as f:
    lines = f.readlines()
    print('total lines:',len(lines))

loop_N = ( len(lines) // (int(atomNum)+2) )
print('iteration number:',loop_N)

if do_atom_colored:
    print('#######  colored begin  #######')
    atom_colored()
    print('#######  colored end  #######\n\n')

if do_atom_to_center_time:
    time.sleep(1)
    print('#######  finding center begin  #######')
    atom_to_center_time()
    print('#######  finding center end  #######\n\n')

if do_center_metal_distribution:
    time.sleep(1)
    print('#######  P(center-metal) begin  #######')
    center_metal_distribution()
    print('#######  finding center end  #######')


