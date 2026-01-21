#!/bin/bash
match_dir=$1
dir=$2
#过滤掉噬菌体
touch $dir/potential_virus_contigs_1.tsv 
touch $dir/potential_virus_contigs_2.tsv 

awk -F"\t" '{if ($16 !~ /;/ && $17 ~ /Viruses/) print}' $match_dir/matches.tsv > $dir/potential_virus_contigs_temp.tsv
# grep -v ";" $match_dir/matches.tsv | grep -i Viruses  > $dir/potential_virus_contigs_temp.tsv
if [ -s $dir/potential_virus_contigs_temp.tsv ]; then
    cut -f16 $dir/potential_virus_contigs_temp.tsv |taxonkit reformat -I 1 -r Unassigned -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}"| sed '1i\TaxID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tStrain'|sed '1d' > $dir/taxa.txt
    paste $dir/potential_virus_contigs_temp.tsv $dir/taxa.txt |cut -f1-3,15,20-  > $dir/potential_virus_contigs_1.tsv

    cat $dir/potential_virus_contigs_1.tsv > $dir/potential_virus_contigs_all.tsv
    rm $dir/potential_virus_contigs_temp.tsv
    rm $dir/taxa.txt
else
    echo "No Virus matches found in potential_virus_contigs_1.tsv." > $dir/potential_virus_contigs_log.tsv
fi

awk -F"\t" '{if ($16 ~ /;/ && $17 ~ /Viruses/) print}' $match_dir/matches.tsv > $dir/multiple_taxids.tsv
# grep ";" $match_dir/matches.tsv | grep -i Viruses  > $dir/multiple_taxids.tsv
if [ -s $dir/multiple_taxids.tsv ]; then
    echo "multiple_taxids.tsv has matches."
    python scripts/python/contigs/split_multiple_taxids.py $dir
    cut -f16 $dir/split_multiple_taxids.tsv |taxonkit reformat -I 1 -r Unassigned -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}"| sed '1i\TaxID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tStrain'|sed '1d' > $dir/taxa.txt
    paste $dir/split_multiple_taxids.tsv $dir/taxa.txt |cut -f1-3,15,20-  > $dir/multiple_taxids_taxa.txt
    python scripts/python/contigs/select_potential_virus.py $dir

    if [ -s $dir/potential_virus_contigs_1.tsv ]; then
        cat $dir/potential_virus_contigs_1.tsv $dir/potential_virus_contigs_2.tsv > $dir/potential_virus_contigs_all.tsv
    else
        cat $dir/potential_virus_contigs_2.tsv > $dir/potential_virus_contigs_all.tsv
    fi
fi
