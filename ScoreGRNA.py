#!/home/clarkg/anaconda2/ python

import sys
import numpy as np



class OffTargetScore(object):

	def pairwise(self,_positions):
		distances=[]
		mean=0
		for p in range(len(_positions)):
			for l in range(p+1,len(_positions)):
				distances.append(abs(_positions[l]-_positions[p]))
		if distances:
			mean=np.mean(distances)
		return mean

	def MIT(self,gRNA,offtarget):
		weights=[0,0,.014,0,0,.395,.317,0,.389,.079,.445,.508,.613,.851,.732,.828,.615,.804,.685,.583]
		scored=0
		_msites=[]
		gRNA=gRNA.upper()
		offtarget=offtarget.upper()
		for k in range(0,len(weights)):
			if gRNA[k] != offtarget[k]:	
				if scored == 0:
					scored=(1-weights[k])
				else:
					scored=scored*(1-weights[k])
				_msites.append(k+1)
		if _msites:
			mdist=self.pairwise(_msites)
		else:
			mdist=1
		_termTwo=1/float( ((19-mdist)/float(19))*4 + 1)
		if len(_msites):
			_termThree=1/float(len(_msites)**2)
			final=scored*_termTwo*_termThree*100
			return final
		else:
			return 0
	



if __name__ == '__main__':
	#Testing here
	
	gRNA="TATGGTTTCAAGTCCCAGTG"#TGG"
	off_o="TATGCTTTCAAGTCCCAGAG"#AGG"
	off_t="TGTGGTTAGAAGTCCCAGTG"#AGG"

	inst=OffTargetScore()

	fnl=inst.MIT(gRNA,off_o)
	print fnl
	#inst.MIT(gRNA,"TGTGGTTAGAAGTCCCAGTGAGG")		
	#inst.MIT(gRNA,"TAAGGTATTGAGTCCCAGTGAGG")
