#!/bin/bash
dir=$1
non_virus=$2
endogenous=$3

cut -f23 $dir/matches.tsv |taxonkit reformat -I 1 -r Unassigned -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}"| sed '1i\TaxID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tStrain'|sed '1d' > $dir/taxa.txt
paste $dir/matches.tsv $dir/taxa.txt |cut -f1-5,20,24- > $dir/matches_taxa.tsv

#grep -vi "virus" $dir/matches_taxa.tsv > $non_virus
#grep -iE "retro|endogenous" $dir/matches_taxa.tsv > $endogenous

non_virus_results=$(grep -vi "virus" $dir/matches_taxa.tsv)
endogenous_results=$(grep -iE "retro|endogenous" $dir/matches_taxa.tsv)

if [ -n "$non_virus_results" ]; then
	  echo "$non_virus_results" > $non_virus
  else
	    echo "No non-virus matches found."
fi

if [ -n "$endogenous_results" ]; then
	  echo "$endogenous_results" > $endogenous
  else
	    echo "No endogenous matches found."
fi
