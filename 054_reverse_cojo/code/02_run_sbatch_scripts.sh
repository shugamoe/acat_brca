#|-- jobs/clumped_cojo
#|-- jobs/plink_clumping
#`-- jobs/raw_cojo

echo "Submitting plink clumping jobs"
for job in jobs/plink_clumping/*; do
  sleep .5s
  echo ${job}
  sbatch ${job}
done
echo "plink clumping jobs submitted"

echo "Submitting raw cojo jobs (no clumping index SNPs)"
for job in jobs/raw_cojo/*; do
  sleep .5s
  echo ${job}
  sbatch ${job}
done
echo "Raw cojo jobs submitted"

echo "Submitting clumped cojo jobs (dependent on outputs of plink clumping)"
for job in jobs/clumped_cojo/*; do
  sleep .5s
  echo ${job}
  sbatch ${job}
done
echo "clumped cojo jobs submitted"
