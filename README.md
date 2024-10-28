<div align="left">
<img src="./docs/logo.png" width="300" height="auto">
</div>

# oncotreeVIS – An interactive graphical user interface  for visualising mutation tree cohorts

<br/>
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

        <table><tr><td>
        The expected input is a JSON file with the following key values:<p style='margin:6px;'></p>
            <table width=650px>
            <tr style='border-bottom: 1px solid darkgray; border-top: 1px solid darkgray'>
            <td><b>Key</b></td><td style='padding-left:20px'><b>Data structure</b></td></tr>
            <tr><td style='vertical-align:top'><b>trees</b></td>
            <td style='padding-left:20px; text-align:justify; vertical-align:top'>
                Nested data structure representing a tree in JSON format, as used in D3.js and anytree (python) libraries. 
                Each node has one or more child nodes (<i>node.children</i>), except for the leaves. In addition, each node 
                has the following attributes: <i>node_id</i> (string/int, required), </i>matching_label</i> (required), 
                <i>size_percent</i> (float, optional), <i>metadata</i> (dictionary, optional), <i>gene_events</i> (dictionary, optional), 
                <i>is_neutral</i> (boolean, optional). The gene_events attribute has two predefined keys (\"mutation\" and \"CNA\"), 
                but any other key names can be used. The values for the \"CNA\" event key are specifically interpreted as amplification 
                or deletion amounts w.r.t. the neutral states. The first three letters of the event key are used 
                in the visualization for displaying a summary for the gene events.<p style='margin:6px;'></p>
                Example of JSON tree: \"AML-03-001\": {\"tree\": {\"node_id\": 0, \"matching_label\": 0, \"children\": 
                [{\"node_id\": 407, \"matching_label\": 14, \"size_percent\": 0.228, \"gene_events\": {\"FLT3-ITD\": {\"mutation\": \"\"}}, 
                \"children\": [{\"node_id\": 408, \"matching_label\": 5, \"size_percent\": 0.772, \"gene_events\": {\"NPM1\": {\"mutation\": 
                \"p.L287fs\"}}}]}]}, \"metadata\": {\"Chemo\": \"No\", \"Gender\": \"Female\", \"VitalStatus\": \"Dead NOS\", \"age\": 59, 
                \"Diagnosis\": \"AML\", \"Response\": \"CR\"}} <p style='margin:6px;'></p>
                Examples for gene_events dictionary: \"gene_events\": {\"NPM1\": {\"mutation\": \"p.L287fs\"}, 
                \"AKT3\": {\"CNA\": 2}, \"JAK2\": {\"CNA\": -1}, \"TP53\": {\"CNA\”: \"-\", \"expression\": \"0.34\"}}
                <p style='margin:6px;'></p>
                </td></tr>
            <tr><td style='vertical-align:top'><b>clusters</b></td>
                <td style='padding-left:20px'>List of lists of tree ids (strings).<p style='margin:6px;'></p>
                Example: [['AML-55-001', 'AML-33-001', 'AML-57-001', 'AML-11-001'],  ['AML-77-001'], ['AML-50-001', 'AML-102-001'], … ]
                <p style='margin:6px;'></p>
                </td></tr>
            <tr><td style='vertical-align:top'><b>pairwise_distances</b></td><td style='padding-left:20px'>
                List of dictionaries where the keys are the tree ids of the pair of trees (strings) and the distance score (float).<p style='margin:6px;'></p>
                Example: [{'sample_1': 'AML-73-001', 'sample_2': 'AML-22-001', 'distance': 0.6072}, … ]
                <p style='margin:6px;'></p>
                </td></tr>
            <tr style='border-bottom: 1px solid darkgray;'><td style='vertical-align:top'><b>highlighted_genes</b></td><td style='padding-left:20px'>
                Styles used: bold, italic, uppercase, color code.<p style='margin:6px;'></p>
                Example: {\"bold\": [\"JAK2\", \"PTEN\", \"TP53\", \"AKT1\"]}
                <p style='margin:6px;'></p>
                </td></tr>
          </table>
          <br/>The JSON files used for the predefined datasets are available on 
          <a href='https://github.com/cbg-ethz/oncotreeVIS/tree/main/data' target=git>github</a>.
        </td><td style='padding:50px'><b>Load JSON file:</b><p style='margin:6px;'></p><input type=file id='picker' onchange='load()'</td></tr></table>

