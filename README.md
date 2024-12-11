# viral-benchmarking
This is a temporary home for the viral assembly benchmarking script.

# Dependencies
- Nextflow
- Podman, Docker, Singularity, or Apptainer

# Test
Below are instructions for how to test the benchmarking tool. The test includes a reference-based assembly generated using the Twist Synthetic Influenza H1N1 (2009) RNA control (read more [here](https://www.twistbioscience.com/products/ngs/synthetic-viral-controls)). The assemblies used to create these synthetic controls were used as "truth" sequences (NC_026431, NC_026432, NC_026433, NC_026434, NC_026435, NC_026436, NC_026437, and NC_026438). The test sample was not particularly high quality, thus variable results are expected.
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
Your results should look like the table below:
|sample                            |nc_errors | slen| tlen| sub| del| ins| miss| other| term| completeness|  accuracy| n_reads| per_reads_mapped| read_cov| mean_read_depth| mean_read_q| mean_map_q|
|:---------------------------------|:---------|----:|----:|---:|---:|---:|----:|-----:|----:|------------:|---------:|-------:|----------------:|--------:|---------------:|-----------:|----------:|
|H1N1-1_T1_Influenza_A-seg-7_MP-6  |null      |  964|  982|   0|  25|   7|    0|     0|   25|     97.45418|  96.65622|   18922|                9|  97.4542|             114|        36.0|       60.0|
|H1N1-1_T1_Influenza_A-seg-5_NP-4  |null      | 1507| 1497|   1|   0|  10|    1|     0|    0|     99.93320|  99.26471|   18922|               12| 100.0000|             112|        35.9|       60.0|
|H1N1-1_T1_Influenza_A-seg-6_NA-2  |null      | 1381| 1410|   0|  37|   8|   73|     0|   37|     92.19858|  96.53846|   18922|                4|  97.3759|              38|        36.3|       60.0|
|H1N1-1_T1_Influenza_A-seg-4_HA-9  |null      | 1591| 1701|   2|   2|   0|   84|     0|  110|     88.59494|  99.73457|   18922|               12|  93.5332|              96|        36.0|       60.0|
|H1N1-1_T1_Influenza_A-seg-8_NS-3  |null      |  791|  863|   0|  14|   0|  103|     0|   72|     79.72190|  97.96512|   18922|                3|  94.6698|              48|        36.1|       59.9|
|H1N1-1_T1_Influenza_A-seg-2_PB1-2 |null      | 2054| 2274|   0| 107|   0|   94|     0|  220|     86.19173|  94.54082|   18922|               11|  91.2929|              64|        36.1|       60.0|
|H1N1-1_T1_Influenza_A-seg-3_PA-1  |null      | 2105| 2151|   0|   0|   0|  127|     0|   46|     91.95723| 100.00000|   18922|               14|  97.8615|              92|        36.1|       60.0|
|H1N1-1_T1_Influenza_A-seg-1_PB2-3 |null      | 2230| 2280|   0|   0|   0|  158|     0|   50|     90.87719| 100.00000|   18922|               34|  98.2456|             201|        36.0|       60.0|
