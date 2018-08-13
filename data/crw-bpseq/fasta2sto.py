import sys

with open(sys.argv[1]) as f:
	lines = f.readlines()
	output = """# STOCKHOLM 1.0
	
{0}{1}
{2}{3}
//	
	"""
	print output.format('template'.ljust(25), lines[1].strip(), '#=GC SS_cons'.ljust(25), lines[2].strip())
	