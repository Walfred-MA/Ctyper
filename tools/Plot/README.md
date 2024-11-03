# README for typemutant.py

This README provides an overview of how to use `typemutant.py` and describes its output figure in detail.

### Python Dependencies

1. **Python 3.7+**
2. **NumPy**
3. **Pandas**  
   [https://pandas.pydata.org/](https://pandas.pydata.org/)
4. **Matplotlib**  
   [https://matplotlib.org/stable/install/index.html](https://matplotlib.org/stable/install/index.html)


## Parameters

Supported parameters:

- **Help:**
  - `-h`: Print help information.
- **Inputs:**
  - `-i <string>`: Path to an individual input file.
- **samples highlight:**
  - `-n <string>`: name of samples highlight.
- **Gff annotation:**
  - `-g <string>`: Path to gff3 annotation file.


### Usages walkthrough
Visualization is performed on a **gene-by-gene** basis (not genome-wide).

For example, to visualize the gene **AMY1A**:


1. **Extract the annotation for AMY_group1 from the full annotation table:**

   ```bash
   grep "^AMY_group1_" PangenomeAlleles_annotationfix.tsv > AMY_group1_annotationfix.tsv
   ```

2. **Extract the genotyping results for AMY_group1:**

   ```bash
   grep "^result: AMY_group1_" genotype.txt
   ```

   Output:

   ```
   result: AMY_group1_GW00031_h1_556,AMY_group1_GW00051_h2_891,
   ```

   The genotyping result is: `AMY_group1_GW00031_h1_556,AMY_group1_GW00051_h2_891`.

3. **Visualize the results:**

   ```bash
   python typemutant.py -i AMY_group1_annotationfix.tsv -n "AMY_group1_GW00031_h1_556,AMY_group1_GW00051_h2_891," -o output.png
   ```

**Optional:** To visualize the GENCODE annotation on the MSA:

4. **Obtain the GENCODE annotation:**

   ```bash
   grep "gene_name=AMY" genecode.gff3 > AMY.gff3
   ```

5. **Run the visualization with GENCODE annotation:**

   ```bash
   python typemutant.py -i AMY_group1_annotationfix.tsv -g AMY.gff3 -n "AMY_group1_GW00031_h1_556,AMY_group1_GW00051_h2_891," -o output.png
   ```

### figure decription

This output will create a mutant map for pengnome-alleles involve.  
Each row on this mutant map corresponds to each pengnome-allele as well as variants on it.   
The colors correspond to different allele-types.   
The red dots on the leftside indicate the location highlighted alleles on the map. 


