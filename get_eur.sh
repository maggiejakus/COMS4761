for i in {1..22};do
 plink2 --pfile all_phase3_ns --chr $i --keep eur_ids --max-alleles 2 --snps-only just-acgt --maf 0.001 --make-bed --out eur.$i
done

