author=Malik

for i in {1..22};do
  echo "CHR BP RSID A1 A2 SE BETA P ESS" | tr ' ' '\t' > ~/Documents/CGNMX/GWAS/${author}/${author}.${i}
  #split on CHROMOSOME
  cat ~/Documents/CGNMX/GWAS/${author}/${author}.txt | awk -v var="$i" '$1== var {print $0}'>> ~/Documents/CGNMX/GWAS/${author}/${author}.${i}
done


