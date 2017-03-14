import sys,math
import numpy as np



def main():
	atomos = PDB2atoms(sys.argv[1])
	ntermA, n_pointer = N_termini(atomos)
	ctermA, c_pointer = C_termini(atomos)


	if len(ntermA) != len(ctermA):
		print "\n\n\n\n\n\t\tERROR !!! there is a different number of \"C\" and \"N\" terminals"
		
	newPDB,n = [],0	
		

	for i in atomos:
		n+=1
		flag = True
		for j in xrange(len(n_pointer)):
			if i==n_pointer[j]:
				flag = False
				for k in ntermA[j]:
					#newPDB.append(k)
					k.fix_atomNum(n)
					k.print_atom()
					n+=1
				#newPDB.append(i)	
				i.fix_atomNum(n)
				i.print_atom()

			if i==c_pointer[j]:
				flag = False
				#newPDB.append(i)
				i.fix_atomNum(n)
				i.print_atom()
				n+=1
				for k in ctermA[j]:
					#newPDB.append(k)	
					k.fix_atomNum(n)
					k.print_atom()

		if flag:
			#newPDB.append(i)	
			i.fix_atomNum(n)
			i.print_atom()


		#n=1
		#for i in newPDB:
			#i.print_atom()
			#j = atom(n, i.symbol, i.residue, i.chain,i.resId, i.atomType, (i.x, i.y, i.z) )
			#j.print_atom()	
			#n+=1


class atom(object):
    
	def __init__(self, atomNum, atomSymbol, residue, chain,resId,atomType, coord): 
		self.symbol = atomSymbol
		self.residue = residue
		self.chain = chain
		self.resId = int(resId)
		self.x = coord[0]
		self.y = coord[1]
		self.z = coord[2]
		self.coord = np.matrix((self.x,self.y,self.z))
		self.atomType = atomType
		self.atomNum = atomNum

	def __eq__(self, other):
		if self.atomNum == other.atomNum and self.symbol == other.symbol and self.residue == other.residue and self.chain == other.chain\
		and self.x == other.x and self.y == other.y and self.z == other.z and self.atomType == other.atomType:
			return True
		else:
			return False

	def print_atom(self):
		print '{:6s}{:5d} {:^4s} {:3s} {:s} {:3d}{:12.3f}{:8.3f}{:8.3f}'.format('ATOM  ',int(self.atomNum),self.symbol,self.residue,self.chain,self.resId,self.x,self.y,self.z)


	def axis(self,other):
		return self.coord - other.coord

	def fix_atomNum(self,newNumber):
		self.atomNum = newNumber



def PDB2atoms(covetz):
	atoms = []
	
	for line in open(covetz):

		# avoid multiple models	
		if line[:5] == 'MODEL':
			if int(line[12:14]) == 1:
				SINGLE_MODEL = True
			else:
				SINGLE_MODEL = False		

		if line[0:6] == 'ATOM  ' and SINGLE_MODEL:
			
			#skip: when two rotamers appear, such as AVAL and BVAL we ignore the BVAL
			if line[16:17]!=' ' and line[16:17]!='A':
				continue	
    		#skip on multimodel. see 1EAW Res 60 chain A
			if line[26:27]!=' ':
				continue

			atomNum, atomSymbol, residue, chain = line[6:11].strip(), line[12:16].strip(),line[17:20].strip(),line[21:22].strip()	
			resId, x, y, z, atomType = line[22:26].strip(),line[30:38],line[38:46],line[46:54],line[77:78].strip()	

			#Inicializo la matriz de objetos atomos				
			atoms.append( atom(atomNum, atomSymbol, residue, chain,resId, atomType, (float(x), float(y), float(z))) )

	return atoms



def rotation(v, axis, theta):

    # Applying the Euler-Rodrigues formula. Angle in degrees.
    axis = np.asarray(axis)
    theta = np.asarray(theta*math.pi/180)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    matrix = np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                       [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                       [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

    return np.dot(matrix, v)

    '''example of use

	v = [3, 5, 0]
	axis = [4, 4, 1]
	theta = 1.2 

	print rotation(v,axis,theta) 

	'''




def N_termini(atoms):
	chain = 'NOTHINGYET'
	ACE, N_pointer = [],[]

	for i in atoms:
		if i.chain != chain:
			N_resNum = i.resId
			chain = i.chain

		if i.chain == chain and N_resNum == i.resId and i.symbol == '3H':
			ACE.append( N_capping_group(i) )
			#NterminalAtoms.append( atom(i.atomNum, i.symbol, i.residue, i.chain,i.resId, i.atomType, (i.x, i.y, i.z)) )

		if i.chain == chain and N_resNum == i.resId and i.symbol == 'N':	
			N_pointer.append(atom(i.atomNum, i.symbol, i.residue, i.chain,i.resId, i.atomType, (i.x, i.y, i.z)) )

	return ACE, N_pointer		


def C_termini(atoms):
	NME, C_pointer = [],[]

	for i in xrange(len(atoms)-1):
		if atoms[i].chain != atoms[i+1].chain or i == len(atoms)-2:
			n = i
			while atoms[n].chain == atoms[i].chain:
				if atoms[n].symbol == 'OXT':
					NME.append( C_capping_group(atoms[n]) )
					#CterminalAtoms.append( atom(atoms[n].atomNum, atoms[n].symbol, atoms[n].residue, atoms[n].chain,atoms[n].resId, atoms[n].atomType, (atoms[n].x, atoms[n].y, atoms[n].z)) )

				if atoms[n].symbol == '2HE2':
					C_pointer.append(atom(atoms[n].atomNum, atoms[n].symbol, atoms[n].residue, atoms[n].chain,atoms[n].resId, atoms[n].atomType, (atoms[n].x, atoms[n].y, atoms[n].z)) )
				n-=1

	return NME, C_pointer

def N_capping_group(N_atom):
	ACE = []
	ACE.append( atom(int(N_atom.atomNum)-3, 'C', 'ACE', N_atom.chain,int(N_atom.resId)-1, 'C', (N_atom.x, N_atom.y, N_atom.z)) )	
	ACE.append( atom(int(N_atom.atomNum)-2, 'O', 'ACE', N_atom.chain,int(N_atom.resId)-1, 'O', (N_atom.x+0.573, N_atom.y-1.064, N_atom.z+0.353)) )	
	ACE.append( atom(int(N_atom.atomNum)-1, 'CH3', 'ACE', N_atom.chain,int(N_atom.resId)-1, 'C', (N_atom.x-1.405, N_atom.y+0.022, N_atom.z-0.499)) )
	return ACE

def C_capping_group(C_atom):
	NME = []
	NME.append( atom(int(C_atom.atomNum)+1, 'N', 'NME', C_atom.chain,int(C_atom.resId)+1, 'N', (C_atom.x, C_atom.y, C_atom.z)) )	
	NME.append( atom(int(C_atom.atomNum)+2, 'CH3', 'NME', C_atom.chain,int(C_atom.resId)+1, 'C', (C_atom.x+1.307, C_atom.y+0.62, C_atom.z+0.328)) )	
	return NME




if __name__ == '__main__':
  main()