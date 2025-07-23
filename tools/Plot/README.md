# README for typemutant.py

This README provides an overview of how to use `typemutant.py` and describes its output figure in detail.

### Python Dependencies

1. **Python 3.7+**
2. **Pandas**
3. **NumPy**
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
  - `-G <string>`: Path to gff3 annotation file.
- **Gff annotation:**
  - `-g <string>`: the gene or gene group name.

### Usages walkthrough
Visualization is performed on a **gene-by-gene** basis (not genome-wide).

For example, to visualize the gene **SMN**:


1. **Visualize the results:**

   ```bash
   python typemutant.py -i PangenomeAlleles_annotationfix.tsv -g SMN -n genotype.txt -o output.png
   ```

**Optional:** To visualize the GENCODE annotation on the MSA:

4. **Obtain the GENCODE annotation:**

   ```bash
   grep "gene_name=SMN" genecode.gff3 > SMN.gff3
   ```

5. **Run the visualization with GENCODE annotation:**

   ```bash
   python typemutant.py -i PangenomeAlleles_annotationfix.tsv -g SMN -G SMN.gff3 -n genotype.txt -o output.png
   ```

### figure decription

1. This output will create a mutant map for pengnome-alleles involve.  
2. Each row on this mutant map corresponds to each pengnome-allele as well as variants on it.   
3. The colors correspond to different allele-types.   
4. The red dots on the leftside indicate the location highlighted alleles on the map.  




