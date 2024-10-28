<div align="left">
<img src="./docs/logo.png" width="300" height="auto">
</div>

# oncotreeVIS – An interactive graphical user interface  for visualising mutation tree cohorts

<br/>

### Web app
----------
<b>A visualization demo can be seen here:</b> https://cbg-ethz.github.io/oncotreeVIS


### Abstract
----------

<p align="justify">In recent years, developments in next-generation sequencing technology and computational methodology have made it possible to reconstruct, with increasing precision, the evolutionary history of tumors and their cell phylogenies, represented as mutation trees. Many mutation tree inference tools exist, yet they only provide visualizations of individual output trees, with limited amount of details, which makes it difficult to inspect the mutation trees at the cohort level and understand their differences and commonalities, an important task in computational oncology highly relevant for tumor board clinical decisions.</p>

<p align="justify">We introduce <b>oncotreeVIS</b>, an interactive graphical user interface for visualizing mutation tree cohorts and tree posterior distributions obtained from mutation tree inference tools. OncotreeVIS can display mutation trees that encode single or joint genetic events, such as point mutations and copy number changes, and highlight subclones with matching mutation events, conserved trajectories, drug-gene interactions, and k-nearest neighbor trees. OncotreeVIS facilitates the visual inspection of mutation tree clusters and pairwise tree distances, if provided by the user. It is available both as a JavaScript library that can be used locally or as a web application that can be accessed online to visualize seven default datasets of public mutation tree cohorts, or new user data provided in a predefined JSON format.</p>
<br/>
<div align="center">
<img src="./docs/poster.png" width="90%" height="auto">
</div>

### Input format
----------
The expected input is a JSON file with the following key values:<br/><br/>

| Key      | Data structure |
| ----------- | ----------- |
| clusters|List of lists of tree ids (strings).<br/><br/>Example:[['AML-55-001', 'AML-33-001', 'AML-57-001', 'AML-11-001'],  ['AML-77-001'], ['AML-50-001', 'AML-102-001'], … ] |
| pairwise_distances   | List of dictionaries where the keys are the tree ids of the pair of trees (strings) and the distance score (float).<br/><br/> Example: [{'sample_1': 'AML-73-001', 'sample_2': 'AML-22-001', 'distance': 0.6072}, … ]|

<br/><br/>The JSON files used for the predefined datasets are available in the <i>data</i> folder. 

