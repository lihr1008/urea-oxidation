# ChemExtractor

## -1 Dependencies
**Requires PYTHON=3.7 and JAVA**

```bash
pip install -r requirements.txt
```

## 0. Prepare Data and Configure Settings

### 0.1 Prepare Data
Currently, only data downloaded from the Web of Science (WOS) website is supported.
Place the downloaded files in the directory [root]/[origin_papers_dn], where [root] and [origin_papers_dn] are defined in the configuration settings (see Section 0.2).

### 0.2 onfigure Settings
Modify the config.json file to set up your project directory. Update the commented values under "proj_config" and leave others unchanged
```json
{
  "proj_config": {
    "name": "# Fill in a project name (optional)",
    "root": "# Fill in the path to the data storage folder",
    "target_mol_type": "# Fill in the target chemical type, e.g., '\''ALLOY'\'' or '\''ORG'\''",
    "struc_info": "ELEMENTS",
    "keywords": "# Fill in keywords (optional)",
    "needed_rows": "# If '\''target_mol_type'\'' is '\''ALLOY'\'', fill in required periodic table rows (format: list, e.g., [3,4,5])",
    "origin_papers_dn": "wos",
    "papers_fn": "papers.tsv",
    "texts_fn": "texts.tsv",
    "texts_tagged_fn": "texts_tagged.tsv",
    "tag_token_pairs_fn": "tag_token_pairs.tsv",
    "mols_fn": "mols.tsv",
    "unique_mols_fn": "unique_mols.tsv",
    "org_mols_fn": "org_mols.tsv",
    "org_mols_img_fn": "org_mols_img.xlsx",
    "elements_count_fn": "elements_count.tsv",
    "elements_link_count_fn": "elements_link_count.tsv",
    "paper_and_mol_fn": "paper_and_mol.tsv",
    "words_vec_fn": "words_vec.txt",
    "cner_vec_fn": "cner_vec.tsv",
    "cner_cluster_fig_fn": "cner_cluster.pdf",
    "cner_cluster_fn": "cner_cluster.tsv",
    "fig_dn": "fig",
    "log_fn": "log.txt"
  }
}
```

## 1. Data Cleaning
Clean and organize data from different sources into a standard format.
```bash
# For PDF files, run:
python data_cleaner/pdf_paper_extractor.py
# For WOS TSV files, run:
python data_cleaner/paper_extractor.py

# After processing either PDF or WOS files, run:
python data_cleaner/text_extractor.py
```

## 2. Text Preprocessing
Tokenization, POS tagging, and NER:
```bash
java -jar rxn-extractor-restful-0.0.3-SNAPSHOT.jar
python text_preprocess/nlp_tagger.py # --> texts_tagged.tsv
python text_preprocess/tag_token_extractor.py # --> tag_token_pairs.tsv
# python text_preprocess/mol_extractor.py
```

## 3. Molecule Extraction
```bash
python mol_extractor/keywords_extractor.py # --> words_vec.txt & keywords_synonyms.tsv
python mol_extractor/mol_extractor.py # --> unique_mols.tsv
python mol_extractor/mol_filter.py # --> filter_mols.tsv

# For alloys, run:
python mol_extractor/alloy_extractor.py # --> elements_count.tsv & elements_link_count.tsv
# For organic compounds, run:
python mol_extractor/org_extractor
```
