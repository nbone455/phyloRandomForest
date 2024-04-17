# phyloRandomForest
This R script was a part of my Macroevolution graduate level class. It looks at how Random forest algorithms can be impacted by phylogenetic traits. It's designed to help people understand the impact of tree-structured trait correlation on classification/ The script investigates the impact of removing different types of traits or tree information on the accuracy of the classifier.

## Model Preparation



### Random Forest Models
Full Model: Utilizes all available features.
No Random Trait Model: Excludes an arbitrary trait.
No Clade Info Model: Excludes clade information.
No Tree-Correlated Trait Model: Excludes a tree-correlated trait.

## Error Analysis
OOB Error: Out-of-bag error is calculated for each model.
Overall Error: Generalized error rates are calculated for test data.

## Visualization
Box plots are used for visualizing the error rates across different models.

Feel free to clone, modify, and use this code for your research and applications.
