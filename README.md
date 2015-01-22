### Text
- [ ] Abstract
  - [ ] Background
  - [ ] Findings
  - [ ] Conclusions
- [ ] Data description
  - [ ] Background
  - [ ] Assembly QC
    - [ ] Percentage of mapped reads
    - [ ] Run ALE or LAP or ?
  - [x] KEGG stuff
    - [x] From gene to KO to Module/Pathway
    - [ ] Explain and refer to Figure 1
  - [ ] Coverage analysis
    - [ ] Explain and refer to Figure 2
  - [ ] Do we need a taxonomic profile?
    - [ ] Cite some old paper instead, working on the same sample
    - [ ] Alternatively, taxator-tk profile, weighted contigs
- [ ] Discussion
- [ ] Update GitHub URL

### Figures
- [x] Prepare Figure 1: metagenome/metatranscriptome
  - [x] http://www.kegg.jp/kegg-bin/show_pathway?org_name=map&mapno=00680&mapscale=&show_description=show
  - [x] http://www.kegg.jp/kegg-bin/show_module?M00567
  - [x] http://www.kegg.jp/kegg-bin/show_module?M00356
  - [x] http://www.kegg.jp/kegg-bin/show_module?M00357
  - [x] http://www.kegg.jp/dbget-bin/www_bget?ec:2.8.4.1
- [ ] Discussion with Alex and Andreas

### Tables
- [x] Update Table 2
- [x] Update Table 3

### Makefile
- [ ] Include KEGG analyses
- [ ] Merge everything into a single makefile
  - [ ] Create different tasks for assembly, mapping, annotation, ...
- [ ] Final run (after updating binaries)
- [ ] Remove CeBiTec-specific stuff (e.g. explicit qsub commands)
- [ ] Eventually commit it to GitHub
