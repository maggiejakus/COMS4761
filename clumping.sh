#rm keep.ss

#this script performs clumping using plink
#clumps the data from GWAS of interest to remove SNPs that are i nhigh LD:

#for each chromosome


#get all SNPs identified in michailidou study, save in summary statistics (ss)
#zcat ~/athena/doc_score/raw_ss/Michailidou/chr_ss/michailidou_${i}.ss.gz > ss
#zcat ~/athena/doc_score/raw_ss/Christophersen/chr_ss/christophersen_${i}.ss.gz > ss

author=christophersen
upper_author=Christophersen

for i in {1..22};do
  cat ~/Documents/CGNMX/GWAS/${upper_author}/${upper_author}.${i} > ss 

  #finds duplicates
  #cat ~/athena/refs/1000genomes/eur.${i}.bim | cut -f2 | sort | uniq -d > dup_ids
  cat ~/Documents/BMEN/eur_files/eur.${i}.bim | cut -f2 | sort | uniq -d > dup_ids

  #perform clumping with plink. exclude duplicates. use 1000 genomes as reference bfile. clump based on rsID
  plink --seed 1 --memory 4000 --threads 1 --bfile ~/Documents/BMEN/eur_files/eur.${i} --exclude dup_ids --clump ss --clump-snp-field RSID --clump-field P --clump-p1 0.1 --clump-r2 0.5 --out temp_out.normal

  #cleans up temp_out.clumped, giving the good SNPs only
  sed -e 's/ [ ]*/\t/g' temp_out.normal.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 > done_rsids

  cat ss | fgrep -w -f done_rsids  > keep.ss.${i}

  #super clump
  plink --seed 1 --memory 4000 --threads 1 --bfile ~/Documents/BMEN/eur_files/eur.${i} --exclude dup_ids --clump ss --clump-snp-field RSID --clump-field P --clump-p1 0.000001 --clump-r2 0.25 --out temp_out.super

  #cleans up temp_out.clumped, giving the good SNPs only
  sed -e 's/ [ ]*/\t/g' temp_out.super.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 > done_rsids

  #fgrep -w -f done_rsids ss >> superclump.${author}.ss.${i} 
  cat ss | fgrep -w -f done_rsids > superclump.${author}.ss.${i}

done

#rm superclump.snps.${author}.to_plot
for i in {1..22};do
  cat superclump.${author}.ss.${i} | cut -f3 >> superclump.snps.${author}.to_plot #CHANGE HERE-cut to variant_id
done
