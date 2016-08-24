# GeneD3

[![Click for interactive visualization with GeneD3](https://github.com/esteinig/geneD3/blob/master/title_image.png)](http://htmlpreview.github.io/?https://github.com/esteinig/geneD3/master/title_image.html)

Interactive sequence comparisons in [D3](https://github.com/d3), orchestrated with Python.

Click the image above to see an interactive example of GeneD3:
* zoom and drag
* hover over blocks to see annotation (GenBank)
* click and hold comparisons to see alignment and identity (BLAST)

### Usage
---

`geneD3.py --files Sequence1.gbk,Sequence2.gbk --names Seq1,Seq2`

### Output
---

The program outputs the visualization as HTML, which you can open in Chrome or Firefox.

You can use [SVG Crowbar](http://nytimes.github.io/svg-crowbar/) to download the visualization from your browser as a static vector file - works in both Chrome and Firefox!

### Dependencies
---

* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [Python 3](https://www.continuum.io/downloads)

* [Pandas](http://pandas.pydata.org/)
* [BioPython](http://biopython.org/wiki/Documentation)

### Installation
---

```
git clone https://github.com/esteinig/geneD3.git $HOME/geneD3
chmod +x $HOME/geneD3/geneD3.py
echo "PATH=$PATH:$HOME/geneD3" >> $HOME/.bashrc
source $HOME/.bashrc
```

### Options
---

`geneD3.py --help`

#### Basics

```
--files           Comma-delimited list of annotation and sequence files, currently only 
                  formatted as complete GenBank
--names           Comma-delimited list of unique names for each sequence
--features        Comma-delimited list of features to parse from annotation (CDS)
--qualifiers      Comma-delimited list of qualifiers to parse from annotations (gene,product)
--path            Project path (.)
--project         Name of project (GeneD3)
--write_segments  Write segment files to project path for manual manipulation,
                  e.g. colours or annotation tooltips (False)
--read_segments   Read segment files from project path, overrides previous
                  properties of segments (False)
```

#### BLAST

```
--blast_path      Path to blastn or blastp executable (blastn)
--makedb_path     Path to makeblastdb (makeblastdb)
--db_type         BLAST DB type, corresponding to executable (nucl)
--identity        Minimum identity to include comparison (0.80)
--length          Minimum alignment length to include comparison (1000)
--evalue          BLAST e-value (1e-06)
```

#### Visualization

```
--title                 Plot title (GeneD3)
--width                 Canvas width (1800)
--height                Canvas height (1200)
--scale                 Scale segments, default to width of canvas (1800)

--gap                   Gap between sequence segments (300)
--block_height          Height of segment blocks (50)
--block_colour          Base colour of segment blocks (darkgray)
--line_width            Width of sequence line (5)
--line_colour           Colour of sequence line (gray)
--polygon_gap           Gap between comparison polygon and segments (10)
--polygon_fill          Polygon colour (darkgray)
--polygon_opacity       Polygon opacity (0.8)
--polygon_line_colour   Polygon line colour (darkgray)
--polygon_line_width    Polygon line width (1)
```


#### Text

```
--title_colour          Title colour (#404040)
--title_opacity         Title opacity (0.8)
--title_font_family     Title font family (times)
--title_font_weight     Title font weight (bold)
--title_font_size       Title font size (400%)

--name_colour           Name colour (#404040)
--name_opacity          Name opacity (0.8)
--name_font_family      Name font family (times)
--name_font_weight      Name font weight (bold)
--name_font_size        Name font size (400%)
```

##### Coordinates

```
--y_start   Start of segments, Y (250)
--x_start   Start of segments, X (500)
--y_text    Segments names, relative to top of segments, Y (10)
--x_text    Segment names, X (0
--y_title   Title, Y (100)
--x_title   Title, defaults to half of canvas width X (900)
```

### Update
---

In the next version:
* pip install on PyPI
* tests and annotations
* download accessions from NCBI
* forward / reverse marker style
* alignment of segments on canvas
* comparison opacity scaled by identity
* additional annotation and sequence formats

