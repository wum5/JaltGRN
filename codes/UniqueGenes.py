import sys


if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: python UniqueGenes.py infile outfile"
		sys.exit()
		
	infile = open(sys.argv[1], "r")
	outfile = open(sys.argv[2], "w")
        geneList = []
	
	for line in infile:
            line = line.rstrip()
            if 'Solyc' not in line: continue
            geneID = (line.split(',')[1]).split('_')[1]
            geneID = geneID.rstrip('"')
            geneList.append(geneID)

        bf = len(geneList)
        geneList = list(set(geneList))
        af = len(geneList)
        ratio = float(af)/bf
        print "Total Genes:%d, Uniq Genes:%d, percent:%f" % (bf, af, ratio)
        
        for key in geneList:
            outfile.write(key+'\n')
