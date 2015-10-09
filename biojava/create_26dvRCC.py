import RCCpackage.RCCobject as rcco
import sys


if __name__ == '__main__':

	PDB = sys.argv[1]
	CHAIN = sys.argv[2]
	if len(sys.argv) >= 4:
		DIR = sys.argv[3]
		rccs = rcco.RCC(PDB,CHAIN,tmpdir=DIR)
	else:
		rccs = rcco.RCC(PDB,CHAIN)

	print ','.join(map(str,rccs.RCCvector))



