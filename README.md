# cpDNA survey
Analyze NCBI cpDNA data for relation between species closeness and chloroplast genome sequence lengths.

Pipeline is used in the research:
"Survey of chloroplast sequence data: how close is close?".

Pipeline is implemented as a stand-alone Python3 script.
Additional script `survey_figures.py` creates images, used in the article, from the data and results.


## Description

Pipeline acquires needed data and stores all the research results.
The tool is controlled by arguments that determine the scope of data to be analyzed and the thresholds used in the calculations.
The only mandatory argument is taxa name, others have default values.
The results of an analysis are saved in an Excel spreadsheet and related figures can be generated.

Workflow:
* Download RefSeq sequence summaries of given taxa. Summary is the same as list that can be generated on
  [NCBI Genome page](https://www.ncbi.nlm.nih.gov/genome/browse#!/organelles/).
* Calculate sequence length distributions for families and genera. Box plot method is used, with calculated properties:
  the median of the sample, the first and third quartiles, and the lower and upper whiskers, and the outliers.
  As the main descriptor for the width of the distribution, ratio of IQR and median values (IQR/median) is used, presented as percentages.
* Examination of the wide distributions. Thresholds for distribution to be declared wide are set with an argument
* Checking of detected outlier sequences.
  For species of outlier sequence, script downloads summary of all species complete chloroplast sequences from
  [NCBI GenBank database](https://www.ncbi.nlm.nih.gov/genbank/).
  Script checks are lengths of acquired sequences closer to the median of genus distribution than the original RefSeq sequence.


## Output

Pipeline stores files in the current working directory. Prefix of filenames is determined by argument values.
List of stored files with used suffixes:
* Excel file, suffix `results.xlsx`, with data of each step stored in one or more worksheets.
* Text summary file, suffix `summary.txt`, stores a lot of unstructured calculated data.
* Figure data, suffix `figures_data.json`, stores data that is used for figure creation.
* Figures with suffixes `data.png`, `iqr_median.png`, `wide_families.png`, and `wide_genera.png`.

In addition to stored files, pipeline uses a directory for caching data downloaded from NCBI.
Directory name is set with an argument. Json format is used for data caching.


## Usage

```
# Help
python3 survey.py -h

# Analyze sequences of one or more taxa
python3 survey.py <taxa> [<taxa>]*

# Command used for the research
python3 survey.py asterids rosids -c asterids -c rosids -p 2022-10-22

# Figure creation
# With main script, add -i switch
python3 survey.py asterids rosids -c asterids -c rosids -p 2022-10-22 -i
# With survey_figures.py script. Figure descriptions correspond to Figures in the article.
python3 survey_figures.py <figures_data.json filename> {data|im|wf|wg}
```

## Requirements

Required libraries:
* [BioPython](https://biopython.org/)
* [numpy](https://numpy.org/)
* [pandas](https://pandas.pydata.org/)
* [scipy](https://scipy.org/)
* [ETE3](http://etetoolkit.org/)
* [Matplotlib](https://matplotlib.org/) (for figures, not required for the analysis)


## Citation

Turudić, A.; Liber, Z.; Grdiša, M.; Jakše, J.; Varga, F.; Šatović, Z.
Variation in Chloroplast Genome Size: Biological Phenomena and Technological Artifacts.
Plants 2023, 12, 254. https://doi.org/10.3390/plants12020254
