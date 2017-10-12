ScoreGRNA.py is a python simple class showing the use of the 'MIT Score' as indicated  from http://crispr.mit.edu/about

For any paired gRNA sequence and off-target, this class returns an individual score serving as a proxy for specificity.

To obtain a general score for the gRNA, we loop over all individual target scores. 
If we assume a [list] of off-targets (case-aware of 'grnaSequence', ACTG-alphabet,without PAM sequence),
the final aggregate score is 'mitScore'

		"""
		from ScoreGRNA import OffTargetScore

		calcS=OffTargetScore()

		def MITcrispr(foundOffs,grnaSequence):
			tally=0
			for offtarg in foundOffs:
				tally+=calcS.MIT(grnaSequence,offtarg)
			mitScore=100/(100+tally)*100
			return mitScore
		"""
