#!usr/bin/perl

### WARNING: remove A/T and G/C (or T/A and C/G) SNPs before running this on your data!

open( OUTBIM, ">", "bimout.bim" );
open( OUTVCF, ">", "vcfout.vcf" );
open( PROBLEM,">", "problemsnps");

$c1=0;$c2=0;$c3=0;$c4=0;$c5=0;$total=0;
while( <> )
{
	$total++;
	chomp;
	@la = split( /	/,$_ );
	$a1r = $la[4];
	$a2r = $la[5];
	$a1r =~ tr/ACGT/TGCA/;
	$a2r =~ tr/ACGT/TGCA/;
	if( $la[4] eq $la[6] )		## Ref is same as 1st allele: just print
	{
		$c1++;
		print OUTBIM $la[0]."\t".$la[1]."\t".$la[2]."\t".$la[3]."\t".$la[4]."\t".$la[5]."\n";
		print OUTVCF "chr".$la[0]."\t".$la[3]."\t".$la[1]."\t".$la[4]."\t".$la[5]."\t.\t.\t.\n";
	}
	elsif( $la[5] eq $la[6] )	## Ref is same as 2nd allele: switch and print
	{
		$c2++;
		print OUTBIM $la[0]."\t".$la[1]."\t".$la[2]."\t".$la[3]."\t".$la[4]."\t".$la[5]."\n";
		print OUTVCF "chr".$la[0]."\t".$la[3]."\t".$la[1]."\t".$la[5]."\t".$la[4]."\t.\t.\t.\n";
	}
	elsif( $a1r eq $la[6] )	## Ref is rev-comp of 1st allele: flip and print
	{
		$c3++;
		print OUTBIM $la[0]."\t".$la[1]."\t".$la[2]."\t".$la[3]."\t".$a1r."\t".$a2r."\n";
		print OUTVCF "chr".$la[0]."\t".$la[3]."\t".$la[1]."\t".$a1r."\t".$a2r."\t.\t.\t.\n";
	}
	elsif( $a2r eq $la[6] )	## Ref is rev-comp of 2nd allele: flip, switch and print
	{
		$c4++;
		print OUTBIM $la[0]."\t".$la[1]."\t".$la[2]."\t".$la[3]."\t".$a1r."\t".$a2r."\n";
		print OUTVCF "chr".$la[0]."\t".$la[3]."\t".$la[1]."\t".$a2r."\t".$a1r."\t.\t.\t.\n";
	}
	else
	{
		$c5++;
		print PROBLEM $la[1]."\n";
	}
}

print STDERR "All done! Totals:\nRef matched allele1 :\t$c1\nRef matched allele2 :\t$c2\nRef matched comp of allele1 :\t$c3\nRef matched comp of allele2 :\t$c4\nProblem SNPs :\t$c5\nTotal SNPs analysed :\t$total\n";

close( OUTBIM );
close( OUTVCF );
close( PROBLEM );
