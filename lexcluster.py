#!/usr/bin/env python3.9
"""
Finding all examples in the Blust PPH reflex list which have:
- A limited geographical distribution;
and/or
- Limited subgroup distribution
August 2022

Zorc's axes (All quotes from Zorc 2021):
    - North Bisayan axis: "[...] West Bisayan, Central Bisayan, Asi’, Bikol, and Hanunoo
(South Mangyan) subgroups as well as Kagayanen. Geographically this includes northern
Panay, southern Mindoro, Tablas, Caluya, Sibale, Banton, Masbate, the Bikol peninsula, and
Cagayancillo Island."
    - Southern Luzon axis: "This axis affects Tagalog and languages in the Central Luzon
subgroup which include Sambalic and Kapampangan), Remontado/Hatang-Kayi, Casiguran Dumagat
Agta, Bulalakawnon, Bikol, and Hanunoo."
    - Eastern Mindanao axis: "This represents a sequence of unique distributions among
South Bisayan, Mamanwa, Mansakan, Danao, Subanen, and various members of Manobo, mostly
located in the eastern area of Mindanao (except Subanen which clearly had moved west). It
was this axis that led Gallman (1997) to propose a “Proto East Mindanao” (which is
demonstrably an axis rather than a subgroup)."
    - Northern Luzon axis: "This represents innovative intrusions among Ilokano, Cagayan
Valley, Central Cordilleran, and Bashiic languages."
    - Central Luzon axis: "This represents cognates mainly drawn from the Central and
South Cordilleran subgroups, with encroachments further south into Sambalic, Kapampangan,
and Tagalog."
    - Palawan-Kalamian axis: "This is another axis limited to an insular area, Palawan and
the Calamian Islands, with the Kalamianic and Palawanic subgroups. Kalamianic is an
isolate, but Palawanic may be part of the Greater Central Philippine subgroup proposed by
Blust (1991). The co-existence of these two very different language groups in the greater
Palawan area appears to have led to these local innovations. Data for the Palawanic
languages comes from Thiessen (1980)."
    - Palawan-Mindoro axis: "This axis formed as a result of relations between the peoples
of Palawan and Mindoro. It involves the North and South Mangyan, Kalamianic, and Palawanic
subgroups."
    - Southern Mindanao axis: "This axis involves interchanges among the Southwestern
Mindanao languages (Tboli, Koranadal and Sarangani Bilaan, Tiruray ~ Teduray),
Southeastern Mindanao (Bagobo Klata), and members of the Greater Central Philippine
subgroup such as Tagalog (when pre-Tagalog was still this far south)57, Danao, and
Manobo."
    - Central Mindanao axis: "This axis involves intimate connections among the Mindanao
members of the Greater Central Philippine subgroup such as Mansakan, Danao, Subanen, and
Manobo, but which did not cross over into Bisayan, Bikol, or Tagalog."
    - Bisayan-Bikol axis: "The outreach of this axis includes most Bisayan and Bikol
languages. Normally, innovations appearing in Bisayan and Bikol languages would be
assigned to their immediate proto-language (Proto Central Philippine). Instead, an axis is
proposed when it is clear that their sister languages (Tagalic and Mansakan) also have
innovations within the same semantic parameters. Note that there is an informal coinage or
portmanteau of “Bisaya” and “Bikol”—Bisakol—which is usually the term for the three
Bisayan languages spoken within the Bicol Region [the Central Bisayan speech varieties of
Central Sorsoganon (Sorsogon), Masbate (Masbatenyo), and Southern Sorsoganon (Gubat)]
which have fallen under the influence of either Northern or Southern Bikol. There is even
one Austronesian reconstruction where cognates seem to be limited to just this single axis
within the Philippines."
    - Catanduanes axis: "The island of Catanduanes (depicted in the map immediately above
and enlarged in the map below) offers the requisite features for an axis relationship. The
northern part has Northern Catanduanes (McFarland’s “Pandan”) which is an isolate (a
subgroup all by itself) while Viracnon (Southern Catanduanes) is a member of the Northern
Bikol subgroup (with Naga and Legaspi across the strait). Despite the genetic distance
between these languages, six shared innovations have come to light as a result of
investigating a 400-word list developed by 26 McFarland during his dissertation research
in 1972. While there is a possibility that these were borrowed by Viracnon from Northern
Catanduanes, the innovations have clearly remained on the island and there is no direct
evidence of any dominant relationship. Surely, the replacement of basic vocabulary such as
‘tooth’ and ‘go home’ is significant"
    - Western Luzon axis: "This axis became reasonably clear to me when various entries in
McFarland’s magnificent survey of Northern Philippine languages (1977) were primarily on
the left side of any given page."
    - North-Central Luzon axis: "Similar to the above, this axis became apparent when
McFarland’s survey of Northern Philippine languages (1977) had cognate entries primarily
on the right side of any given page.  Careful examination of the map below illustrates
that these do not actually cover “Eastern”, but rather “North-Central” areas of
Luzon. These are the only four that have thus far survived scrutiny."
    - Bashiic-Ilokano axis: "In my review (Zorc 2020) of Blust (2019), I rejected
out-of-hand 19 of the reconstructions by Blust (2019) that only had Bashiic (1) and
Ilokano (2a) cognates as evidence for a PPH etymon. All of these can be found in the
ACD. We do not know the prehistory of the region, but careful consideration of the
cognates, reflexes, and distribution indicates that such forms may have been innovated in
the northernmost region of the Philippines, but did not take hold and carry on over or
down into the Cordilleran and more southerly Philippine languages. I find these
reconstructions analogous to the situation of Tagalog páwis and Kapampangan páwas for
‘sweat’ (discussed above in 2.2.). While several of these have the “right” and
differential reflexes (such as for *R), these etyma should be treated as evidence for an
axis relationship until other cognates are discovered to raise these to PNP or PPH."
    - The last axis is defined phonologically: the *R > Y language axis / the North
extension.
----> 13 clusters
"""
import csv
import itertools
import sys
import tabulate

import geopy.distance as distance
import shapely.geometry as geometry

from cldfcatalog import Config
from pyglottolog import Glottolog
from pickle import dump, load
from pathlib import Path

MICROGROUPS = [
    "Batanic",
    "Northern Luzon",
    "Central Luzon",
    "Ati",
    "Kalamian",
    "Bilic",
    "South Mangyan",
    "Palawanic",
    "Central Philippine",
    "Manobo",
    "Danaw",
    "Subanen",
    "Sangiric",
    "Minahasan",
    "Gorontalo-Mongondow"
]

MGCODES = [
    "bata1315", # Batanic
    "nort3238", # Northern Luzon
    "nort2873", # Northern Mangyan
    "sang1335", # Sangiric
    "bili1253", # Bilic
    "cent2080", # Central Luzon
    "kala1389", # Kalamian
    "mina1272", # Minahasan
    "cent2080", # Central Luzon
    "atii1237", # Ati
    "lamp1241", # Lampungic
    "mano1276", # Manobo
    "pala1354", # Palawanic
    "goro1257", # Gorontalo/Mongondow
    "cent2246", # Central Philippine
    "suba1253", # Subanen
    "dana1253", # Danaw
    "umir1236", # Umiray Dumaget Agta
]

## The 7 most spoken languages in the Philippines according to the 2010 Philippines
## census, and major regional lingua francae cf. Ascuncion-Lande (1971)
REGIONALS = [
    "taga1270", # Tagalog
    "ilok1237", # Ilokano
    "pamp1243", # Kapampangan
    "biko1240", # Bikol (Language family. need to find out what the "prestige" language is
    "cebu1242", # Cebuano / Bisayas
    "hili1240", # Hiligaynon
    "wara1300" # Waray
]
    
class GlottoCache:
    """Save loaded Glottolog languoids so we don't have to wait every time"""
    def __init__(self, path, fname="glottocache.pickle"):
        if type(path) == str:
            path = Path(path)
        self.cachepath = path.joinpath(fname)
        try:
            with open(self.cachepath, "rb") as cachefile:
                self.cache = load(cachefile)
        except FileNotFoundError:
            self.cache = {}

        cfg = Config.from_file()
        self.glottolog = Glottolog(cfg.get_clone("glottolog"))
        # So we know whether to save
        self.written = False

    def __del__(self):
        self.save()

    def get(self, glottocode):
        if not glottocode:
            return None
        elif glottocode in self.cache:
            return self.cache[glottocode]
        else:
            lg = self.glottolog.languoid(glottocode)
            self.cache[glottocode] = lg
            self.written = True
            return lg

    def save(self):
        if self.written:
            with open(self.cachepath, "wb+") as cachefile:
                dump(self.cache, cachefile)


def load_data(path, delimiter="\t"):
    with open(path) as f:
        rows = [row for row in csv.DictReader(f, delimiter=delimiter)]
    return rows


def attach_glottolog_data(glottocache, rows):
    """Add location and subgroup data from Glottolog"""
    cache = {}
    result = []
    for row in rows:
        code = row["GlottoCode"]
        if code: # Some ACD entries don't have glottocodes
            if code not in cache:
                lg = glottocache.get(code)
                cache[code] = {"Latitude": lg.latitude,
                               "Longitude": lg.longitude,
                               "Ancestors": [a.name for a in lg.ancestors]}
            result.append(row | cache[code])
    return result


def interpolate_positions(glottocache, rows):
    """Estimate missing coordinate data.
    Some languages as referred to in the ACD are listed as families in Glottolog,
    e.g. Bikol. Treat the positions for these as the centroid of the positions
    of the member languages of the family
    """
    cache = {}
    result = []
    for row in rows:
        code = row["GlottoCode"]
        latitude = row["Latitude"]
        lg = glottocache.get(code)
        if latitude == None and lg.category == "Family":
            if code not in cache:
                lg = glottocache.get(code)
                member_coords = [(m.latitude, m.longitude) for m in lg.iter_descendants() if m.latitude]
                centroid = geometry.Polygon(member_coords).centroid
                cache[code] = centroid
            else:
                centroid = cache[code]
            try:
                row["Latitude"] = centroid.x
                row["Longitude"] = centroid.y
                row["InterpolatedDistance"] = True    
            except IndexError:
                print(row["GlottologName"], centroid.wkt)
        else:
            row["InterpolatedDistance"] = False
        result.append(row)
    return result


def groupby(rows, field):
    """Group rows according to field value. Returns {field_value: rows} dict"""
    grouped = {}
    for row in rows:
        datum = row[field]
        if datum in grouped:
            grouped[datum].append(row)
        else:
            grouped[datum] = [row]
    return grouped


def compute_distances(rows):
    """Calculate summary statistics for distances between languages which have
    reflexes of this cognate set
    """
    coords = [(row["Latitude"], row["Longitude"]) for row in rows] # geopy.distance.distance wants latitude first
    # Distance between each pair of languages
    distances = []
    for a, b in itertools.product(coords, repeat=2):
        try:
            if a != b:
                km = distance.distance(a, b).km
                distances.append(km)
        except ValueError:
            print("Couldn't calculate distance for {}, {}".format(a, b))
    return set(distances)


def get_unique_microgroups(glottocache, groups, rows, attr="glottocode"):
    """Count unique subgroups for this cognate set"""
    microgroups = []
    for row in rows:
        lg = glottocache.get(row["GlottoCode"])
        ancestors = [getattr(a, attr) for a in lg.ancestors]
        for group in groups:
            if group in ancestors:
                microgroups.append(group)
    return set(microgroups)


def summarise_lexical_data(glottocache, grouped):
    """Search for subgroup or geographically limited cognate sets.
    Instead of using an arbitrary distance cutoff, calculate the maximum distance
    between reflexes within a set.
    """
    result = []
    for protoform, rows in grouped.items():
        if len(rows) > 1:
            # Distances
            distances = compute_distances(rows)
            unique_groups = get_unique_microgroups(glottocache, MICROGROUPS, rows, attr="name")
            set_row = {
                "protoform": protoform,
                "reflexes": len(rows),
                "maxdist": max(distances),
                "mindist": min(distances),
                "meandist": sum(distances) / len(distances),
                "interpolated": True if any([row["InterpolatedDistance"] for row in rows]) else False,
                "microgroups": unique_groups,
                "nmicrogroups": len(unique_groups),
                "hasregionallang": has_languages(REGIONALS, rows)
            }
            result.append(set_row)
    return result


def has_languages(languages, rows):
    """Utility function to add a column for languages of interest e.g. linguas franca"""
    return bool(len([r for r in rows if r["GlottoCode"] in languages]))

##
## Clustering analysis: attempting to confirm Zorc's (2021) interaction axes by
## applying a clustering algorithm to the ACD Philippines data. Try several different
## algorithms from scikit-learn and see which gives the best results
## September 2022
##

import numpy
import pandas
import sklearn.cluster as skcluster
import gower
import hdbscan

import geopandas
import matplotlib.pyplot as plot
import matplotlib.patches as patches

def build_microgroup_matrix(rows):
    """Transform the result of find_suspicious_sets into a presence/absence
    matrix with one column per microgroup
    """
    matrix = []
    for row in rows:
        matrix_row = {"protoform": row["protoform"], "meandist": row["meandist"]}
        for mg in MICROGROUPS:
            state = 1 if mg in row["microgroups"] else 0
            matrix_row[mg] = state
        matrix.append(matrix_row)
    return matrix


def build_language_matrix(glottocache, rows):
    """| Language | Microgroup | (Lng) | (Lat) | Set_A | ... | Set_N |"""
    matrix = []
    all_sets = set([row["ProtoForm"] for row in rows])
    lang_groups = groupby(rows, "GlottoCode")
    for code, rows in lang_groups.items():
        if code:
            lg = gc.get(code)
            lg_sets = set([row["ProtoForm"] for row in rows])
            matrix_row = {
                "name": lg.name,
                "glottocode": lg.glottocode,
                "x": lg.longitude,
                "y": lg.latitude
            }
            for set_ in all_sets:
                state = 1 if set_ in lg_sets else 0
                matrix_row[set_] = state
            matrix.append(matrix_row)
    return pandas.DataFrame(matrix)

## Some convenience functions for counting cognate set reflexes per language and per set
## Use these to experiment with how filtering sets / languages with low state counts
## the clustering

def get_n_row_feats(matrix, col_regex):
    """How many features are present for each row in a binary matrix?"""
    row_sums = [(n, c) for n, c in zip(matrix.iloc[:,0], matrix.filter(regex=col_regex).sum(axis=1))]
67    return pandas.DataFrame(row_sums)


def get_feat_counts(matrix, col_regex):
    """How many times is each individual feature represented in the data?"""
    return matrix.filter(regex=col_regex).sum(axis=0).sort_values(ascending=False)


def make_clusters(featmatrix, metric, min_size=2):
    frame = pandas.DataFrame(featmatrix)
    if metric == "gower":
        distance_matrix = gower.gower_matrix(frame)
        clusterer = hdbscan.HDBSCAN(metric="precomputed")
        # .astype(numpy.float64): github.com/scikit-learn-contrib/hdbscan/issues/71
        clusterer.fit(distance_matrix.astype(numpy.float64))
    else:
        clusterer= hdbscan.HDBSCAN(metric=metric, min_cluster_size=min_size)
        clusterer.fit(featmatrix)
    # Add computed clusters to feature matrix
    return clusterer, frame


def display_clusters(matrix):
    """Convenient readout so I can see what languages/sets are in what cluster"""
    name_clust = lm[["name", "cluster"]]
    return name_clust.groupby("cluster")["name"].apply(list)


def make_map(matrix):
    """Plot clusters from matrix with clusters added"""
    # Draw map of the Philippines
    phmap = geopandas.read_file(
        "/Users/isaac_stead/Projects/lexcluster/data/philippines-outline.json"
    )
    fig, axis = plot.subplots(figsize=(12, 6))
    phmap.plot(color="lightgrey", ax=axis)
    # Get relevant data out of matrix
    xs = matrix["x"]
    ys = matrix["y"]
    clusters = matrix["cluster"]
    n_clusters = len(set(clusters))
    names = matrix["name"]
    # Assign colours to cluster indices
    colours = plot.cm.rainbow(np.linspace(0, 1, n_clusters))
    # Plot and label each language point coloured by cluster index
    for x, y, clust, name in zip(xs, ys, clusters, names):
        axis.scatter(x, y, color=colours[clust])
        axis.annotate(name, (x, y), size=4)
    # Create the legend
    handles = [patches.Patch(color=colours[clust], label=clust if clust>=0 else "outliers")
               for clust in set(clusters)]
    axis.legend(handles=handles)
    plot.show()
            

# First run, with HDBSCAN fit against Gower distance matrix: 39 clusters, 394 outliers
# Hypothesis: outliers will be cognate sets which appear in most microgroups, i.e. cognate
# sets genuinely descended from some protolanguage?
# Also, this is with the mean distance measure, try running it without the mean distance
# and seeing how this changes the results, as the mean distance between reflexes is maybe
# epiphenomenal -- clusters should emerge from the data without explicitly including the
# distance...?
# Clustering the cognate sets might also be the wrong approach... instead cluster languages
# directly. The feature matrix could look something like:
# | Language | Microgroup | (Lng) | (Lat) | Set_A | ... | Set_N |
# Using this layout of the data and a distance measure appropriate for binary data such as
# Jaccard, Dice etc seems to give better results: recovering 10 clusters with a min
# cluster size of 2. Now I need to write a quick function to tell me what the clusters are
# and then plot them on a map...

###
### Tying it together:
###

def load_and_process_sets(glottocache, infile):
    """Load Philippine cognate set data and return summary table and matrix"""
    data = load_data(infile)
    rows = attach_glottolog_data(glottocache, data)
    rows = interpolate_positions(glottocache, rows)
    grouped = groupby(rows, "ProtoForm")
    summary = summarise_lexical_data(gc, grouped)
    mm = build_microgroup_matrix(summary)
    lm = build_language_matrix(glottocache, rows)
    return lm, mm, summary

def main():
    gc = GlottoCache(".")
    infile, summaryfile, matrixfile = sys.argv[1:]
    matrix, summary = load_and_process_sets(infile)

    with open(summaryfile, "w") as f:
        writer = csv.DictWriter(f, fieldnames=summary[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)
    with open(matrixfile, "w") as f:
        writer = csv.DictWriter(f, fieldnames=matrix[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(matrix)
    

if __name__ == "__main__":
    main()

