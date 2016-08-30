#!/usr/bin/env python

import os
import json
import pandas
import textwrap
import argparse
import subprocess

from Bio import SeqIO
from html.parser import HTMLParser

def main():

    """

    """

    command_line = CommandLine()

    args = command_line.arg_dict

    project_path = os.path.join(args["path"], args["project"])

    os.makedirs(project_path, exist_ok=True)

    sequences = [Sequence(genbank=args["files"][i], name=args["names"][i], features=args["features"],
                          qualifiers=args["qualifiers"]) for i in range(len(args["files"]))]

    cg = ComparisonGenerator(sequences, path=project_path, blast_path=args["blast_path"],
                             makedb_path=args["makedb_path"], db_type=args["db_type"], type="sequential",
                             identity=args["identity"], min_length=args["min_length"], evalue=args["evalue"])

    comparisons = cg.comparisons

    options = {"style": "block",
               "block_height": args["block_height"],
               "block_colour": args["block_colour"],
               "gap": args["gap"],
               "line_width": args["line_width"],
               "line_colour": args["line_colour"],
               "polygon_fill": args["polygon_fill"],
               "polygon_line_colour": ["polygon_line_colour"],
               "polygon_line_width": args["polygon_line_width"],
               "polygon_opacity": args["polygon_opacity"],
               "polygon_gap": args["polygon_gap"],
               "title_colour": args["title_colour"],
               "name_colour": args["name_colour"],
               "title_font_weight": args["title_font_weight"],
               "title_font_family": args["title_font_family"],
               "title_font_size": args["title_font_size"],
               "title_opacity": args["title_opacity"],
               "name_font_weight": args["name_font_weight"],
               "name_font_family": args["name_font_family"],
               "name_font_size": args["name_font_size"],
               "name_opacity": args["name_opacity"],
               "y_start": args["y_start"],
               "y_text": args["y_text"],
               "y_title": args["y_title"],
               "x_start": args["x_start"],
               "x_text": args["x_text"],
               "x_title": args["x_title"]}

    if args["read_segments"]:
        read_files = [os.path.join(project_path, name + ".csv") for name in args["names"]]
    else:
        read_files = None

    if args["write_segments"]:
        write_files = [os.path.join(project_path, name + ".csv") for name in args["names"]]
    else:
        write_files = None

    viz = GeneD3(sequences, comparisons, title=args["title"], height=args["height"], width=args["width"],
                 project=args["project"], align=args["align"],scale=args["scale"], write_segments=write_files,
                 read_segments=read_files, path=project_path, options=options)

    viz.write()


### Command Line ###


class CommandLine:

    """ Command Line Parser """

    def __init__(self):

        self.parser = argparse.ArgumentParser(description='GeneD3 v.0.1', add_help=True)
        self.set_parser()

        self.args = self.parser.parse_args()
        self.arg_dict = vars(self.args)

    def set_parser(self):

        """Initiate command line parsing module"""

        self.parser.add_argument("--files", dest='files', type=lambda s: [str(item) for item in s.split(',')],
                                 required=True, help="List of files (complete GenBank), comma-delimited"
                                                     "(e.g. File1.gb,File2.gb,File3.gb)")
        self.parser.add_argument("--names", dest='names', type=lambda s: [str(item) for item in s.split(',')],
                                 required=True, help="List of sequence names, comma-delimited (e.g. Seq1,Seq2, Seq2)")
        self.parser.add_argument("--features", dest='features', default="CDS",
                                 type=lambda s: [str(item) for item in s.split(',')],
                                 required=False, help="Features to parse and visualize from annotation [CDS]")
        self.parser.add_argument("--qualifiers", dest='qualifiers', default="gene,product",
                                 type=lambda s: [str(item) for item in s.split(',')], required=False,
                                 help="Qualifiers to parse and visualize from annotation [gene,product]")

        self.parser.add_argument('--read_segments', dest='read_segments', action="store_true", default=False,
                                 required=False, help="Write segment files for modification [False]")
        self.parser.add_argument('--write_segments', dest='write_segments', action="store_true", default=False,
                                 required=False, help="Read segment files for modification [False]")

        self.parser.add_argument('--path', dest='path', default=os.getcwd(), required=False, type=str,
                                 help="Project path [.]")
        self.parser.add_argument('--project', dest='project', default='GeneD3', required=False, type=str,
                                 help="Project name [GeneD3]")

        self.parser.add_argument('--blast_path', dest='blast_path', default='blastn', required=False, type=str,
                                 help="Path to BLAST [blastn]")
        self.parser.add_argument('--makedb_path', dest='makedb_path', default='makeblastdb', required=False, type=str,
                                 help="Path to MAKEBLASTDB [makeblastdb]")
        self.parser.add_argument('--db_type', dest='db_type', default='nucl', required=False, type=str,
                                 help="BLAST DB type [nucl]")

        self.parser.add_argument('--identity', dest='identity', default=0.80, required=False, type=float,
                                 help="Minimum identity to keep sequence alignments [0.80]")
        self.parser.add_argument('--length', dest='min_length', default=1000, required=False, type=float,
                                 help="Minimum length to keep sequence alignment [1000]")
        self.parser.add_argument('--evalue', dest='evalue', default=0.000001, required=False, type=float,
                                 help="BLAST e-value [1e-06]")

        self.parser.add_argument('--title', dest='title', default='GeneD3', required=False, type=str,
                                 help="Title [GeneD3]")
        self.parser.add_argument('--width', dest='width', default=1800, required=False, type=int,
                                 help="Canvas width [1800]")
        self.parser.add_argument('--height', dest='height', default=1200, required=False, type=int,
                                 help="Canvas height [1200]")
        self.parser.add_argument('--scale', dest='scale', default=1800, required=False, type=int,
                                 help="Scale segments to maximum length [1800]")

        self.parser.add_argument('--align', dest="align", default="left", required=False, type=str,
                                 help="Alignment of sequences [left]")
        self.parser.add_argument('--style', dest='style', default='block', required=False, type=str,
                                 help="Marker style [block]")

        self.parser.add_argument('--gap', dest="gap", default=300, required=False, type=int,
                                 help="Gap between segments [300]")
        self.parser.add_argument('--block_height', dest='block_height', default=50, required=False, type=int,
                                 help="Block height [50]")
        self.parser.add_argument('--block_colour', dest='block_colour', default=1800, required=False, type=str,
                                 help="Block colour [darkgray]")

        self.parser.add_argument('--line_width', dest='line_width', default=5, required=False, type=int,
                                 help="Line width [5]")
        self.parser.add_argument('--line_colour', dest='line_colour', default='gray', required=False, type=str,
                                 help="Line colour [gray]")

        self.parser.add_argument('--polygon_gap', dest="polygon_gap", default=10, required=False, type=int,
                                 help="Polygon gap to segments [10]")
        self.parser.add_argument('--polygon_fill', dest='polygon_fill', default='gray', required=False, type=str,
                                 help="Polygon fill [gray]")
        self.parser.add_argument('--polygon_line_colour', dest='polygon_line_colour', default='darkgray', required=False,
                                 type=str, help="Polygon line colour [gray]")
        self.parser.add_argument('--polygon_line_width', dest='polygon_line_width', default=2, required=False, type=int,
                                 help="Polygon line width [2]")
        self.parser.add_argument('--polygon_opacity', dest='polygon_opacity', default=0.8, required=False, type=int,
                                 help="Polygon opacity [0.8]")

        self.parser.add_argument('--title_colour', dest="title_colour", default="#404040", required=False, type=str,
                                 help="Title colour [#404040]")
        self.parser.add_argument('--title_font_family', dest="title_font_family", default="times", required=False, type=str,
                                 help="Title font family [times]")
        self.parser.add_argument('--title_font_weight', dest="title_font_weight", default="bold", required=False, type=str,
                                 help="Title font weight [bold]")
        self.parser.add_argument('--title_font_size', dest="title_font_size", default="400%", required=False, type=str,
                                 help="Title font size [400%%]")
        self.parser.add_argument('--title_opacity', dest="title_opacity", default=1, required=False, type=float,
                                 help="Title opacity [0.8]")

        self.parser.add_argument('--name_colour', dest="name_colour", default="#404040", required=False, type=str,
                                 help="Segment names colour [#404040]")
        self.parser.add_argument('--name_font_family', dest="name_font_family", default="times", required=False, type=str,
                                 help="Name font family [times]")
        self.parser.add_argument('--name_font_weight', dest="name_font_weight", default="bold", required=False, type=str,
                                 help="Name font weight [bold]")
        self.parser.add_argument('--name_font_size', dest="name_font_size", default="200%", required=False, type=str,
                                 help="Name font size [200%%]")
        self.parser.add_argument('--name_opacity', dest="name_opacity", default=1, required=False, type=float,
                                 help="Name opacity [0.8]")

        self.parser.add_argument('--y_start', dest="y_start", default=250, required=False, type=float,
                                 help="Y-coordinate for segment start [250]")
        self.parser.add_argument('--y_text', dest="y_text", default=10, required=False, type=float,
                                 help="Y-coordinate for name text relative to segment top [10]")
        self.parser.add_argument('--y_title', dest="y_title", default=100, required=False, type=float,
                                 help="Y-coordinate for title [100]")

        self.parser.add_argument('--x_start', dest="x_start", default=500, required=False, type=float,
                                 help="X-coordinate for segment start [500]")
        self.parser.add_argument('--x_text', dest="x_text", default=0, required=False, type=float,
                                 help="X-coordinate for name of each segment [0]")
        self.parser.add_argument('--x_title', dest="x_title", default=900, required=False, type=float,
                                 help="X-coordinate for title [900]")

### Helper Classes ###

class Tooltip:

    """Tooltip class (under construction), will hold more options to customize Tooltips for D3"""

    def __init__(self):

        self.text_color = '#565051'
        self.head_color = '#565051'

    def get_qualifiers(self, location, qualifiers, order):

        """ Qualifier Tooltip """

        html = '<strong><span style="color:{head}">Location: </span></strong>'.format(head=self.head_color) +\
               '<span style="color:{text}">{start} - {end}</span><br>'.format(text=self.text_color,
                                                                              start=str(location.start),
                                                                              end=str(location.end))

        if len(qualifiers) > 0:
            for qualifier in order:
                try:
                    q = " ".join([word.capitalize() for word in qualifier.split("_")])
                    html += '<strong>' + '<span style="color:' + self.head_color + '">' + q + ": "\
                            + '</span>' + '</strong>' + '<span style="color:' + self.text_color + '">'\
                            + qualifiers[qualifier] + ' ' + '</span>' + '<br>'

                except KeyError:
                    pass

        return html

    def get_blast(self, data):

        htmls = []

        for index, row in data.iterrows():
            string = '<strong>' + '<span style="color:' + self.head_color + '">' + "Reference" + ": " +\
                     '</span>' + '</strong>' + '<span style="color:' + self.text_color + '">' + row["reference"]\
                     + " (" + str(row["start.1"]) + ' - ' + str(row["end.1"]) + ")" + '</span>' + '<br>' +\
                     '<strong>' + '<span style="color:' + self.head_color + '">' + "Query" + ": " +\
                     '</span>' + '</strong>' + '<span style="color:' + self.text_color + '">' + row["query"]\
                     + " (" + str(row["start.2"]) + ' - ' + str(row["end.2"]) + ")" + '</span>' + '<br>' +\
                     '<strong>' + '<span style="color:' + self.head_color + '">' + "Alignment" + ": " +\
                     '</span>' + '</strong>' + '<span style="color:' + self.text_color + '">' + str(row["length"]) +\
                     ' bp' + '</span>' + '<br>' +\
                     '<strong>' + '<span style="color:' + self.head_color + '">' + "Identity" + ": " +\
                     '</span>' + '</strong>' + '<span style="color:' + self.text_color + '">' +\
                     str(row["identity"]) + " %" + '</span>' + '<br>'

            htmls.append(string)

        return htmls

class Blaster:

    def __init__(self, directory=os.getcwd(), blast_path="blastn", makedb_path="makeblastdb", db_type="nucl"):

        self.directory = directory
        self.blast_path = blast_path
        self.makedb_path = makedb_path

        self.db_type = db_type
        self.split = "-vs-"
        self.data = pandas.DataFrame()

        self.comparison = ""

        self.headers = ["query", "reference", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                        "sstart", "send", "evalue", "bitscore"]

    def _make_db(self, reference):

        self.db = os.path.splitext(reference)[0] + ".db"

        with open(os.devnull, "w") as devnull:
            subprocess.call([self.makedb_path, '-in', os.path.join(self.directory, reference),
                             '-dbtype', self.db_type, '-out', os.path.join(self.directory, self.db)],
                            stdout=devnull)

    def run(self, reference, query, identity, length, evalue):

        self._make_db(reference)
        self.comparison = query + self.split + reference

        cmd = [self.blast_path, '-query',  os.path.join(self.directory, query),
               '-db', os.path.join(self.directory, self.db),
               '-outfmt', '6', '-evalue', str(evalue),
               '-out', os.path.join(self.directory, self.comparison)]

        subprocess.call(cmd)

        self.data = pandas.read_csv(os.path.join(self.directory, self.comparison), header=None, sep="\t")
        self.data.columns = self.headers

        self.data["reference"] = reference
        self.data["query"] = query

        self.data = self.data[(self.data["pident"] >= identity*100) & (self.data["length"] >= length)]


### Comparison Classes ###


class ComparisonGenerator:

    def __init__(self, sequences, path=os.getcwd(), blast_path="blastn", makedb_path="makeblastdb",
                 db_type="nucl", type="sequential", identity=0.80, min_length=1000, evalue=1e-06,
                 internal=True):

        self.path = path

        self.sequences = sequences
        self.data = pandas.DataFrame()

        self.db_type = db_type
        self.blast_path = blast_path
        self.makedb_path = makedb_path

        self.identity = identity
        self.min_length = min_length
        self.evalue = evalue

        self.type = type

        self.blaster = Blaster(directory=self.path, blast_path=self.blast_path, makedb_path=self.makedb_path,
                               db_type=self.db_type)

        self.comparisons = list()

        if internal:
            self.run_blast()

    def _write_fastas(self):

        files = []

        for sequence in self.sequences:

            with open(sequence.genbank) as genbank_file:
                sequences = SeqIO.parse(genbank_file, "genbank")
                with open(os.path.join(self.path, sequence.name + ".fa"), "w") as fasta_file:
                    SeqIO.write(sequences, fasta_file, "fasta")

            files.append(sequence.name + ".fa")

        return files

    def _get_order(self, files, order_dict=None):

        if self.type is "sequential":
            sequential = list()
            for i in range(len(files)):
                try:
                    sequential.append((files[i], files[i+1]))
                except IndexError:
                    pass

            return sequential

    def run_blast(self):

        files = self._write_fastas()
        order = self._get_order(files)

        for reference_file, query_file in order:

            reference = os.path.basename(os.path.splitext(reference_file)[0])
            query = os.path.basename(os.path.splitext(query_file)[0])

            self.blaster.run(reference_file, query_file, identity=self.identity, length=self.min_length,
                             evalue=self.evalue)

            comparison = Comparison(comparison=self.blaster.data, reference=reference, query=query)

            self.comparisons.append(comparison)

class Comparison:

    """
    Class holding data for comparison DataFrame and transformation of different comparison types.
    """

    def __init__(self, comparison, reference, query, type="blast"):

        self.type = type
        self.reference = reference
        self.query = query

        self.comparison = comparison

        self.tooltip = Tooltip()

        self.data = pandas.DataFrame()

        self.columns = ["query", "reference", "identity", "length",  "start.2", "end.2", "start.1", "end.1"]

        self.indices = {"blast": [0, 1, 2, 3, 6, 7, 8, 9]}

        self._transform()

    def _transform(self):

        if self.type is "blast":

            if isinstance(self.comparison, pandas.DataFrame):
                self.data = self.comparison[self.indices[self.type]]
            else:
                self.read(self.comparison)

            self.data.columns = self.columns

            tooltips = self.tooltip.get_blast(self.data)

            self.data = self.data.assign(tooltip=tooltips)

    def read(self, file, sep="\t"):

        self.data = pandas.read_csv(file, columns=self.indices[self.type], sep=sep)


### Sequence Classes ###


class Sequence:

    """
    Sequence object holding file paths, sequence information and annotations.
    Download files and parse from accession number at GenBank.
    """

    def __init__(self, genbank, name="Sequence", features=None, qualifiers=None, accession=None):

        self.name = name
        self.genbank = genbank

        self.length = 0

        self.features = features
        self.qualifiers = qualifiers

        self.accession = accession

        self.tooltip = Tooltip()

        self.data = pandas.DataFrame

        self.parse()

    def parse(self):

        """ Parse sequence and annotations from files """

        with open(self.genbank) as gb_file:
            record = SeqIO.read(gb_file, format="genbank")
            self.length = len(record.seq)

            features = [feature for feature in record.features if feature.type in self.features]
            self._transform(features, record)

    def _transform(self, features, record):

        """ Transforms list of annotation features to data frame for input to visualization module """

        rows = []
        for feature in features:
            sequence = {"id": self.name, "length": len(record.seq), "start": int(feature.location.start),
                         "end": int(feature.location.end), "type": feature.type, "strand": feature.strand}

            qualifiers = self._parse_qualifiers(feature)
            tooltip = self._make_tooltip(feature.location, qualifiers)

            rows += [merge_dicts(sequence, qualifiers, tooltip)]

        self.data = pandas.DataFrame(rows, columns=["id", "length",  "type", "start", "end", "strand"]
                                                    + self.qualifiers + ["tooltip"])

    def _parse_qualifiers(self, feature):

        qualifiers = {}

        for qualifier in self.qualifiers:
            try:
                qualifiers[qualifier] = ("; ".join(feature.qualifiers[qualifier]))
            except KeyError:
                pass

        return qualifiers

    def _make_tooltip(self, location, qualifiers):

        return {"tooltip": self.tooltip.get_qualifiers(location, qualifiers, order=self.qualifiers)}


### GeneD3 and Visualization ###


class GeneD3:

    def __init__(self, sequences, comparisons,  options=None, title="GeneD3", height=10800, width=1200,
                 project="geneD3", align="left", scale=1800, write_segments=None, read_segments=None, path=os.getcwd()):

        self.sequences = sequences
        self.comparisons = comparisons

        self.title = title

        self.width = width
        self.height = height

        self.scale = scale
        self.align = align

        self.write_segments = write_segments
        self.read_segments = read_segments

        self.path = path
        self.file = os.path.join(path, project + ".html")

        self.default_options = {"style": "block", "block_height": 50, "block_colour": "darkgray",
                                "gap": 300, "line_width": 5, "line_colour": "gray",
                                "polygon_fill": "gray", "polygon_line_colour": "darkgray",
                                "polygon_line_width": 1, "polygon_opacity": 0.8,
                                "polygon_gap": 10, "title_colour": "#202020",
                                "name_colour": "#202020", "title_font_weight": "bold", "title_font_family": "times",
                                "title_font_size": "400%", "title_opacity": 1, "name_font_weight": "bold",
                                "name_font_family": "times", "name_font_size": "200%", "name_opacity": 1,
                                "y_start": 250, "y_text": 10, "y_title": 100, "x_start": 500, "x_text": 0,
                                "x_title": self.width/2}

        if options is None:
            self.options = self.default_options
        else:
            self.options = self.default_options
            for key, value in options.items():
                try:
                    self.options[key] = value
                except KeyError:
                    print("Could not find specified option:", key)
                    pass

        self.coordinates = pandas.DataFrame()
        self.data = dict()

        self.transform()

    def _get_coordinates(self):

        y = self.options["y_start"]

        coordinates = dict()

        for sequence in self.sequences:

            y_block_top = y - (self.options["block_height"]/2)
            y_block_bottom = y + (self.options["block_height"]/2)
            y_polygon_top = y_block_bottom + self.options["polygon_gap"]
            y_polygon_bottom = y + self.options["gap"] - self.options["block_height"]/2 - self.options["polygon_gap"]

            coords = {"y_line": y, "y_block_bottom": y_block_bottom, "y_block_top": y_block_top,
                      "y_polygon_top": y_polygon_top, "y_polygon_bottom": y_polygon_bottom,
                      "x_start": self.options["x_start"]}

            coordinates[sequence.name] = coords

            y += self.options["gap"]

        self.coordinates = pandas.DataFrame(coordinates)

    def _get_text_data(self, segments):

        segment_text = [{"text": segment.name,
                         "y": self.coordinates.get_value("y_line", segment.name) + self.options["y_text"],
                         "x": self.options["x_text"],
                         "font_weight": self.options["name_font_weight"],
                         "font_family": self.options["name_font_family"],
                         "colour": self.options["name_colour"],
                         "font_size": self.options["name_font_size"],
                         "opacity": self.options["name_opacity"]} for segment in segments]

        title_text = [{"text": self.title,
                       "y": self.options["y_title"],
                       "x": self.options["x_title"],
                       "font_weight": self.options["title_font_weight"],
                       "font_family": self.options["title_font_family"],
                       "colour": self.options["title_colour"],
                       "font_size": self.options["title_font_size"],
                       "opacity": self.options["title_opacity"]}]

        return {"title": title_text, "names": segment_text}

    def transform(self):

        self._get_coordinates()

        segments = [Segment(sequence, coordinates=self.coordinates, block_height=self.options["block_height"],
                            block_colour=self.options["block_colour"], line_colour=self.options["line_colour"],
                            line_width=self.options["line_width"]) for sequence in self.sequences]

        polygons = [Polygon(comparison, coordinates=self.coordinates, fill=self.options["polygon_fill"],
                            line_colour=self.options["polygon_line_colour"],
                            line_width=self.options["polygon_line_width"],
                            opacity=self.options["polygon_opacity"]) for comparison in self.comparisons]

        if self.read_segments:
            for segment in segments:
                segment.read_blocks(os.path.join(self.path, segment.name + ".csv"))

        if self.write_segments:
            for segment in segments:
                segment.write_blocks(os.path.join(self.path, segment.name + ".csv"))

        text = self._get_text_data(segments)

        self.data = {"blocks": [block for seg in segments for block in seg.block_data],
                     "lines": [line for seg in segments for line in seg.line_data],
                     "polygons": [gon for poly in polygons for gon in poly.data],
                     "text": text}

    def write(self):

        options = {"canvas.width": self.width, "canvas.height": self.height, "segment.scale": self.scale}

        viz = Visualization(data=self.data, options=options)

        viz.write_html(file=self.file)

## Blocks and Lines ##

class Segment:

    """
    Segments represent a series of annotations (blocks) along a sequence (line)

    """

    def __init__(self, sequence=None, coordinates=None, block_height=50, block_colour="gray",
                 line_width=10, line_colour="darkgray"):

        self.data = sequence.data
        self.length = sequence.length
        self.name = sequence.name

        self.coordinates = coordinates

        self.block_height = block_height
        self.block_colour = block_colour

        self.line_width = line_width
        self.line_colour = line_colour

        self.blocks = pandas.DataFrame()
        self.line = pandas.DataFrame()

        self.block_data = dict()
        self.line_data = dict()

        self._transform()

    def write_blocks(self, file):

        column_order = ["name", "x", "y", "width", "height", "colour", "note", "text"]

        blocks = self.blocks.copy()

        blocks["note"] = [strip_tags(tooltip) for tooltip in blocks["text"]]

        blocks = blocks[column_order]

        blocks.to_csv(file)

    def read_blocks(self, file):

        self.blocks = pandas.read_csv(file)
        self.block_data = self.blocks.to_dict(orient="records")

    def _transform(self):

        self.blocks["x"] = self.data["start"].add(self.coordinates.get_value("x_start", self.name))
        self.blocks["y"] = self.coordinates.get_value("y_block_top", self.name)
        self.blocks["width"] = self.data["end"].subtract(self.data["start"])
        self.blocks["height"] = self.block_height
        self.blocks["colour"] = self.block_colour
        self.blocks["text"] = self.data["tooltip"]
        self.blocks["name"] = self.data["id"]

        self.line = pandas.DataFrame({"x1": self.coordinates.get_value("x_start", self.name),
                                      "x2": self.length + self.coordinates.get_value("x_start", self.name),
                                      "y1": self.coordinates.get_value("y_line", self.name),
                                      "y2": self.coordinates.get_value("y_line", self.name),
                                      "length": self.length + self.coordinates.get_value("x_start", self.name),
                                      "colour": self.line_colour, "width": self.line_width},
                                     index=pandas.Series(self.name),
                                     columns=["x1", "x2", "y1", "y2", "width", "length", "colour"])

        self.block_data = self.blocks.to_dict(orient="records")
        self.line_data = self.line.to_dict(orient="records")

## Comparison Polygon ##

class Polygon:

    def __init__(self, comparison,  coordinates, fill="gray", line_colour="darkgray", line_width=1, opacity=0.8):

        self.comparison = comparison
        self.reference = comparison.reference
        self.query = comparison.query

        self.data = list()

        self.coordinates = coordinates

        self.line_colour = line_colour
        self.line_width = line_width
        self.fill = fill
        self.opacity = opacity

        self._transform()

    def _transform(self):

        self.data = [{
            "text": row["tooltip"],
            "poly": [{"x": self.coordinates.get_value("x_start", self.reference) + row["start.1"],
                      "y": self.coordinates.get_value("y_polygon_top", self.reference)},
                     {"x": self.coordinates.get_value("x_start", self.reference) + row["end.1"],
                      "y": self.coordinates.get_value("y_polygon_top", self.reference)},
                     {"x": self.coordinates.get_value("x_start", self.reference) + row["end.2"],
                      "y": self.coordinates.get_value("y_polygon_bottom", self.reference)},
                     {"x": self.coordinates.get_value("x_start", self.reference) + row["start.2"],
                      "y": self.coordinates.get_value("y_polygon_bottom", self.reference)},
                     {"x": self.coordinates.get_value("x_start", self.reference) + row["start.1"],
                      "y": self.coordinates.get_value("y_polygon_top", self.reference)}],
            "colour": self.line_colour,
            "width": self.line_width,
            "opacity": self.opacity,
            "fill": self.fill
        } for index, row in self.comparison.data.iterrows()]

## Script ##

class Visualization:

    """Helper class Visualization, holds script for JS D3. Methods to write replace options from Ring Generator in
    script and write the HTML. Initialize with options dict and data from Ring Generator. """

    def __init__(self, data, options):

        self.options = options
        self.data = data

        self.head = textwrap.dedent("""
                    <!DOCTYPE html>
                    <html lang="en">
                    <style>

                        .blockText {
                            fill: #6B6B6B;
                            font-size: 11px;
                            font-family: Courgette, sans-serif;
                        }

                        div.tooltip {
                            position: absolute;
                            text-align: left;
                            max-width: 300px;
                            padding: 11px;
                            font-size: 11px;
                            font-family: Courgette, sans-serif;
                            background: #FEFCFF;
                            border-radius: 11px;
                            pointer-events: none;
                        }

                        .title {
                            fill: #565051;
                            font-size: 14px;
                            font-family: Courgette, sans-serif;
                        }
                    </style>
                    <head>
                        <meta charset="UTF-8">
                        <title></title>
                    </head>
                    <body>
                    <script src="http://d3js.org/d3.v3.min.js" charset="utf-8"></script>

                        """)

        self.data1 = textwrap.dedent("""
                        <script type="application/json" id="blocks">
                     """)

        self.data2 = textwrap.dedent("""
                        </script>
                        <script type="application/json" id="lines">
                     """)

        self.data3 = textwrap.dedent("""
                        </script>
                        <script type="application/json" id="chords">
                     """)

        self.data4 = textwrap.dedent("""
                        </script>
                        <script type="application/json" id="text">
                     """)

        self.script = textwrap.dedent("""
                        </script>
                        <script>
                        var width = canvas.width;
                        var height = canvas.height;
                        var segment_scale = segment.scale;

                        var blocks = JSON.parse(document.getElementById("blocks").innerHTML);
                        var lines = JSON.parse(document.getElementById("lines").innerHTML);
                        var chords = JSON.parse(document.getElementById("chords").innerHTML);
                        var text = JSON.parse(document.getElementById("text").innerHTML);

                        var scaleX = d3.scale.linear()
                                       .domain([0, d3.max(lines, function(d){ return d.length; })])
                                       .range([0, segment_scale]);


                        // Main Body

                        var chart = d3.select("body")
                                    .append("svg")
                                    .attr("id", "Main")
                                    .attr("width", width)
                                    .attr("height", height)
                                    .call(d3.behavior.zoom().on("zoom", function () {
                                        chart.attr("transform", "translate(" + d3.event.translate + ")" + " scale(" + d3.event.scale + ")")
                                    }))
                                    .append("g");

                        // Tooltip CSS

                        var tooltip = d3.select("body").append("div")
                                                       .attr("class", "tooltip")
                                                       .style("opacity", 0);
                        // Line Shell

                        var lineShell = chart.append("g")
                                              //.attr("transform", "translate(" + width / 2 + "," + height / 2 + ")")
                                              .attr("id", "LineShell");

                        lineShell.selectAll("line")
                                 .data(lines)
                                 .enter()
                                 .append("line")
                                 .attr("x1", function(d, i){ return scaleX(d.x1); })
                                 .attr("y1", function(d, i){ return d.y1; })
                                 .attr("x2", function(d, i){ return scaleX(d.x2); })
                                 .attr("y2", function(d, i){ return d.y2; })
                                 .attr("stroke", function(d, i){ return d.colour; })
                                 .attr("stroke-width", function(d, i){ return d.width; })

                        // Block Shell

                        var blockShell = chart.append("g")
                                              //.attr("transform", "translate(" + width / 2 + "," + height / 2 + ")")
                                              .attr("id", "BlockShell");

                        blockShell.selectAll("rect")
                                  .data(blocks)
                                  .enter()
                                  .append("rect")
                                  .attr("id", function(d,i) { return "block_"+i; })
                                  .attr("x", function(d, i){ return scaleX(d.x); })
                                  .attr("y", function(d, i){ return d.y; })
                                  .attr("width", function(d, i){ return scaleX(d.width) })
                                  .attr("height", function(d, i){ return d.height; })
                                  .style("fill", function(d, i){ return d.colour; })
                                  .on("mouseover", function(d){
                                                                tooltip.transition()
                                                                       .duration(200)
                                                                       .style("opacity", .9);
                                                                tooltip.html(d.text)
                                                                       .style("left", (d3.event.pageX + 20) + "px")
                                                                       .style("top", (d3.event.pageY + 10) + "px");
                                                               })

                                  .on("mouseout", function(d) {
                                                                tooltip.transition()
                                                                       .duration(200)
                                                                       .style("opacity", 0)
                                                               });

                        // Comparison Shell

                        var comparisonShell = chart.append("g")
                                              //.attr("transform", "translate(" + width / 2 + "," + height / 2 + ")")
                                              .attr("id", "ComparisonShell");

                        var lineFunction = d3.svg.line()
                                             .x(function(d) { return scaleX(d.x); })
                                             .y(function(d) { return d.y; })
                                             .interpolate("linear");

                        comparisonShell.selectAll("path")
                                  .data(chords)
                                  .enter()
                                  .append("path")
                                   .attr("d", function(d) { return lineFunction(d.poly); } )
                                   .attr("stroke", function(d) { return d.colour; })
                                   .attr("stroke-width", function(d) { return d.width; })
                                   .attr("fill", function(d) { return d.fill; })
                                   .attr("opacity", 0.3)
                                   .on("mouseover", function(){
                                                        d3.select(this).transition().duration(300).style("opacity", 0.8);
                                                    })
                                   .on("mouseout", function(){
                                                d3.select(this).transition().duration(300).style("opacity",  0.3);
                                            })
                                   .on("mousedown", function(d){
                                            tooltip.transition()
                                                .duration(200)
                                                .style("opacity", .9);
                                            tooltip.html(d.text)
                                                .style("left", (d3.event.pageX + 20) + "px")
                                                .style("top", (d3.event.pageY + 10) + "px");
                                            })
                                    .on("mouseup", function(d) {
                                            tooltip.transition()
                                                .duration(200)
                                                .style("opacity", 0)
                                        })

                        var titleShell = chart.append("g")
                                  .attr("id", "TitleShell");

                        titleShell.selectAll("text")
                              .data(text.title)
                              .enter()
                              .append("text")
                              .attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; })
                              .style("text-anchor", "middle")
                              .style("font-size", function(d) { return d.font_size; })
                              .style("font-weight", function(d) { return d.font_weight; })
                              .style("font-family", function(d) { return d.font_family; })
                              .style("fill", function(d) { return d.colour; })
                              .attr("class", "title")
                              .text(function(d) { return d.text; })
                              .style("opacity", function(d) { return d.opacity; });

                        var nameShell = chart.append("g")
                                  .attr("id", "NameShell");

                        nameShell.selectAll("text")
                              .data(text.names)
                              .enter()
                              .append("text")
                              .attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; })
                              .style("text-anchor", "end")
                              .style("font-size", function(d) { return d.font_size; })
                              .style("font-weight", function(d) { return d.font_weight; })
                              .style("font-family", function(d) { return d.font_family; })
                              .style("fill", function(d) { return d.colour; })
                              .attr("class", "name")
                              .text(function(d) { return d.text; })
                              .style("opacity", function(d) { return d.opacity; });


                        </script>

                        </body>
                        </html>

                     """)

        self._set_script()

    def _set_script(self):

        """Replace placeholder values in script with given options."""

        for placeholder, value in self.options.items():
            self.script = self.script.replace(str(placeholder), str(value))

    def write_html(self, file):

        """Write script to HTML."""

        with open(file, 'w') as outfile:
            outfile.write(self.head)

        with open(file, 'a') as outfile:
            outfile.write(self.data1)
            json.dump(self.data['blocks'], outfile, indent=4, sort_keys=True)
            outfile.write(self.data2)
            json.dump(self.data['lines'], outfile, indent=4, sort_keys=True)
            outfile.write(self.data3)
            json.dump(self.data['polygons'], outfile, indent=4, sort_keys=True)
            outfile.write(self.data4)
            json.dump(self.data['text'], outfile, indent=4, sort_keys=True)
            outfile.write(self.script)


### Global Helpers & Tests ###

def merge_dicts(*dict_args):
    
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.

    http://stackoverflow.com/questions/38987/how-can-i-merge-two-python-dictionaries-in-a-single-expression/
    26853961#26853961

    """

    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

class MLStripper(HTMLParser):

    def __init__(self):

        super().__init__()

        self.reset()
        self.strict = False
        self.convert_charrefs= True
        self.fed = []

    def handle_data(self, d):
        self.fed.append(d)

    def get_data(self):
        return ' '.join(self.fed)

def strip_tags(html):
    s = MLStripper()
    s.feed(html)
    return s.get_data()

main()
