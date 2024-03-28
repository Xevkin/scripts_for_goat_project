import sys

derived_cols = sys.argv[1].split(",")

anc_cols = sys.argv[2].split(",")

out_col = sys.argv[3]

with open(sys.argv[4]) as file:

	for line in file:

		line = line.rstrip("\n")

		if line.startswith("chr"):

			continue

		split_line = line.split()

		derived_genotypes = []

		ancestral_genotypes = []

		anc_genotype = split_line[int(out_col)]

		for derived_column in derived_cols:

			derived_genotypes.append(split_line[int(derived_column)])

		for anc_column in anc_cols:

			ancestral_genotypes.append(split_line[int(anc_column)])

		if all(geno == derived_genotypes[0] for geno in derived_genotypes) and (derived_genotypes[0] != anc_genotype) and all(geno2 == anc_genotype for geno2 in ancestral_genotypes):

			print line
