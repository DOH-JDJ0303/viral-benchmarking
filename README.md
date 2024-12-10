# viral-benchmarking
This is a temporary home for the viral assembly benchmarking script.

# Dependencies
- Nextflow
- Podman, Docker, Singularity, or Apptainer

# Test
## 1. Clone the repo
```
git clone https://github.com/DOH-JDJ0303/viral-benchmarking
```
### 2. Run the test samplesheet
Replace `-profile` with the desired container engine.
```
nextflow run viral-benchmarking/main \
    -profile (podman|docker|singularity|apptainer) \
    --input viral-benchmarking/example.csv \
    --outdir results
```
### 3. Check the results
```
cat results/metrics.csv
```
