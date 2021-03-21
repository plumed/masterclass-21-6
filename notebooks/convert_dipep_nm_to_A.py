# This is a short script to convert the diala_traj_nm input that you 
# can find in the mapping/rt39 regtest in PLUMED into A so that it can 
# be analysed in this exercise.  This is what I used to make the 
# trajectory that we have analysed using thse various different 
# variables.

f = open("diala_traj_nm.xyz", "r" )
for i in range(546) : 
    print( int(f.readline()) )
    box = f.readline().split()
    assert( len(box)==3 )
    print( 10*float(box[0]), 10*float(box[1]), 10*float(box[2]) )
    for j in range(22):
        atom = f.readline().split()
        assert( len(atom)==4 ) 
        print( atom[0], 10*float(atom[1]), 10*float(atom[2]), 10*float(atom[3]) )
f.close()
