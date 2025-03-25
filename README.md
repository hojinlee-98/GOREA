# GOREA

## This tool is for summarizing GOBP and extracting meaningful biological context.

In GeneOntology directory, processed databases using GOBP originated from MSigDB (https://www.gsea-msigdb.org/gsea/msigdb) were contained.

To re-construct the database for specific version of MSigDB, utilize the following simple code. 

```bash
python ./scripts/20250124_go_term_id_mapping_hj.py [json file from MSigdb] [output file name]
```
- `20250124_go_term_id_mapping_hj.py` script is included in scripts directory under GeneOntology.

### Example
```bash
python ./scripts/20250124_go_term_id_mapping_hj.py ./v2024.1.Hs/GOBP/c5.go.bp.v2024.1.Hs.json ./v2024.1.Hs/GOBP/20250124_c5.go.bp.v2024.1.Hs.term.id_mapping_hj.txt
```
