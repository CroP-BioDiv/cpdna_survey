#!/usr/bin/env python3

import os.path
import json
import numpy
import pandas
import csv
import subprocess
from collections import defaultdict
from itertools import chain
from statistics import median
from Bio import Entrez
from ete3 import NCBITaxa
ncbi_taxa = NCBITaxa()


def sheets_to_excel(filename, sheets):
    writer = pandas.ExcelWriter(filename)
    for name, columns, rows in sheets:
        df = pandas.DataFrame(rows, columns=columns)
        df.to_excel(writer, sheet_name=name, index=False, header=True)
    writer.save()


def remove_dupl_names(clade_names, last_clade_names, remove_dup_clades):
    lc = clade_names
    if remove_dup_clades and last_clade_names:
        clade_names = [(a if a != b else None) for a, b in zip(clade_names, last_clade_names)]
    return clade_names, lc


def name_2_taxid(name):
    "Check NCBITaxa.get_fuzzy_name_translation() implementation"
    result = ncbi_taxa.db.execute(f'SELECT taxid FROM species WHERE spname = "{name}" LIMIT 1;')
    if taxid := result.fetchone():
        return int(taxid[0])

    result = ncbi_taxa.db.execute(f'SELECT taxid FROM synonym WHERE spname = "{name}" LIMIT 1;')
    if taxid := result.fetchone():
        return int(taxid[0])

    taxid, _, _ = ncbi_taxa.get_fuzzy_name_translation(name)
    return taxid


def num_species_below(taxid):
    print('Num species for', taxid)
    result = ncbi_taxa.db.execute(
        f"SELECT COUNT(*) FROM species WHERE rank = 'species' AND INSTR(track, ',{taxid},') > 0;")
    return int(result.fetchone()[0])


class GroupByTaxonomy:
    class _Taxon:
        def __init__(self, taxid, parent, rank, name):
            self.taxid = taxid
            self.parent = parent     # _Taxon object or None
            self.rank = rank
            self.name = name
            self.sequence_data = []  # List of dicts with extracted data
            self.children = []       # List of _Taxon objects

        lengths = property(lambda self: [s['length'] for s in self.sequence_data])

        def num_of_rank(self, rank):
            if self.rank == rank:
                return 1
            return sum(c.num_of_rank(rank) for c in self.children)

        def num_of_rank_dict(self, rank, data):
            cn = 1 if (self.rank == rank) else sum(c.num_of_rank_dict(rank, data) for c in self.children)
            data[self.taxid] = cn
            return cn

    # Stores data (list of sequences) for given taxa and group arguments
    def __init__(self, sequence_data, filter_clades=None, ranks=('order', 'family', 'genus')):
        self.ranks = ranks
        self.filter_clades = filter_clades

        # Filter sequences with taxid
        self.without_taxid = [d for d in sequence_data if not d['taxid']]
        sequence_data = [d for d in sequence_data if d['taxid']]
        seq_taxids = [d['taxid'] for d in sequence_data]

        # Find all taxids, species and parent clades
        all_taxids = set(seq_taxids)
        # taxid_2_lineage = dict((t, ncbi_taxa.get_lineage(t)) for t in all_taxids)
        taxid_2_lineage = ncbi_taxa.get_lineage_translator(all_taxids)
        all_taxids.update(chain.from_iterable(taxid_2_lineage.values()))
        taxid_2_rank = ncbi_taxa.get_rank(all_taxids)
        taxid_2_name = ncbi_taxa.get_taxid_translator(all_taxids)

        # Extract only taxid of interest
        group_taxids = set(t for t, r in taxid_2_rank.items() if r in self.ranks)
        if filter_clades:
            filter_clades = set(filter(None, map(name_2_taxid, filter_clades)))
            group_taxids.update(filter_clades)

        # Set sequences
        self.taxa = dict()           # taxid -> _Taxon object
        self.top_taxa = set()
        self.parent_taxids = dict()  # taxid -> filtered linage
        for taxid, seq_data in zip(seq_taxids, sequence_data):
            if not (t2l := taxid_2_lineage.get(taxid)):
                print(f'Warning: no linage for taxid {taxid}, species {taxid_2_name.get(taxid)}!')
                continue
            parent_taxids = [t for t in taxid_2_lineage[taxid] if t in group_taxids]
            if not parent_taxids:
                print(f'Warning: no parent taxids for {taxid}, species {taxid_2_name.get(taxid)}!')
                continue
                # raise ValueError(f'No linage info for taxid {taxid}, species {taxid_2_name.get(taxid)}!')

            if filter_clades and filter_clades.isdisjoint(taxid_2_lineage[taxid]):
                parent_taxids = [0] + parent_taxids  # Add Taxon "others"

            # Add nodes. Note: taxon creation is from higher to lower taxonomy
            self.top_taxa.add(parent_taxids[0])
            self.parent_taxids[taxid] = parent_taxids
            taxon = None  # Used as parent (last) taxon
            for taxid in parent_taxids:
                if taxid == 0:  # "others"
                    taxon = self._add_taxon(taxid, None, 'clade', 'others')
                else:
                    taxon = self._add_taxon(taxid, taxon, taxid_2_rank[taxid], taxid_2_name[taxid])
                taxon.sequence_data.append(seq_data)

    def _add_taxon(self, taxid, parent, rank, name):
        if taxon := self.taxa.get(taxid):
            assert taxon.taxid == taxid, (taxon.taxid, taxid, name)
            assert taxon.parent == parent, (taxon.parent, taxid, name)
        else:
            self.taxa[taxid] = taxon = GroupByTaxonomy._Taxon(taxid, parent, rank, name)
            if parent:
                parent.children.append(taxon)  # Add child
        return taxon

    def num_of_rank(self, rank):
        return sum(1 for t in self.taxa.values() if t.rank == rank)

    def num_of_rank_dict(self, rank):
        data = dict()
        for t in self.taxa.values():
            t.num_of_rank_dict(rank, data)
        return data

    def get_top_taxa(self):
        return sorted((self.taxa[i] for i in self.top_taxa),
                      key=lambda x: x.name if x.name != 'others' else '~')  # Set "others" to the end

    def iterate_leaves(self):
        def _il(taxons):
            if (cs := taxons[-1].children):
                for c in self.sorted_children(taxons[-1]):
                    yield from _il(taxons + [c])
            else:
                yield taxons

        for t in self.get_top_taxa():
            yield from _il([t])

    def _il(self, taxons, rank):
        if taxons[-1].rank == rank:
            yield taxons
        else:
            for c in self.sorted_children(taxons[-1]):
                yield from self._il(taxons + [c], rank)

    def iterate_by_rank(self, rank):
        for t in self.get_top_taxa():
            yield from self._il([t], rank)

    def group_by_rank(self, rank):
        return [(t, [taxons[-1] for taxons in self._il([t], rank)]) for t in self.get_top_taxa()]

    #
    def sorted_children(self, taxon):
        return sorted(taxon.children, key=lambda t: (-self.ranks.index(t.rank), t.name))

    def get_table_data(self, row_data_method, sort_method, remove_dup_clades):
        columns, _clade_names = self.column_data()
        rows = []

        def _rows(taxons, last_clade_names):
            taxon = taxons[-1]
            if taxon.children:
                for child in self.sorted_children(taxon):  # sorted(taxon.children, key=lambda t: t.name):
                    last_clade_names = _rows(taxons + [child], last_clade_names)

            else:
                clade_names, last_clade_names = remove_dupl_names(_clade_names(taxons), last_clade_names, remove_dup_clades)
                for seq_data in sorted(taxon.sequence_data, key=sort_method):
                    rows.append(clade_names + row_data_method(seq_data))
                    if remove_dup_clades:
                        clade_names = [None] * len(clade_names)
            return last_clade_names

        for taxon in self.get_top_taxa():
            _rows([taxon], None)
        return columns, rows

    def get_taxon_data(self, taxon_data_method, remove_dup_clades):
        columns, _clade_names = self.column_data()
        rows = []

        def _rows(taxons, last_clade_names):
            taxon = taxons[-1]
            to_continue, row_data = taxon_data_method(taxon)

            if 'N' not in to_continue:
                clade_names, last_clade_names = remove_dupl_names(_clade_names(taxons), last_clade_names, remove_dup_clades)
                rows.append(clade_names + row_data)

            if ('S' not in to_continue) and taxon.children:
                for child in self.sorted_children(taxon):  # sorted(taxon.children, key=lambda t: t.name):
                    last_clade_names = _rows(taxons + [child], last_clade_names)
            return last_clade_names

        for taxon in self.get_top_taxa():
            _rows([taxon], None)
        return columns, rows

    def column_data(self):
        if self.filter_clades:
            columns = ['Clade'] + [r.capitalize() for r in self.ranks]
            num_c = len(columns)

            def _clade_names(taxons):
                r = [None] * num_c
                for t in taxons:
                    r[(self.ranks.index(t.rank) + 1) if t.rank in self.ranks else 0] = t.name
                return r
        else:
            columns = [r.capitalize() for r in self.ranks]
            num_c = len(columns)

            def _clade_names(taxons):
                r = [None] * num_c
                for t in taxons:
                    r[self.ranks.index(t.rank)] = t.name
                return r
        return columns, _clade_names


class _NumSpeciesBelow:
    def __init__(self, cache_dir):
        self.cache_dir = cache_dir
        self.cache_file = os.path.join(cache_dir, 'species_below.json')
        if os.path.exists(self.cache_file):
            with open(self.cache_file, 'r', encoding='utf-8') as _in:
                self._data = dict((int(k), int(v)) for k, v in json.load(_in).items())
        else:
            self._data = dict()

    def get(self, taxid):
        taxid = int(taxid)
        if not (num := self._data.get(taxid)):
            num = num_species_below(taxid)
            self._data[taxid] = num
        return num

    def save(self):
        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)
        with open(self.cache_file, 'w', encoding='utf-8') as _out:
            json.dump(self._data, _out, indent=2)


class _Results:
    def __init__(self, params):
        self.sequence_data = fetch_refseq(params)  # Fetch RefSeq data
        self.gbt_l = GroupByTaxonomy(self.sequence_data, filter_clades=params.clade,
                                     ranks=('order', 'family', 'subfamily', 'tribe', 'subtribe', 'genus'))
        self.gbt_s = GroupByTaxonomy(self.sequence_data, filter_clades=params.clade,
                                     ranks=('order', 'family', 'genus'))
        self.taxid_2_bp = self.get_taxid_2_bp(self.gbt_s)
        self.taxid_2_num_genera = self.gbt_s.num_of_rank_dict('genus')
        self.genus_taxid_2_outliers = self.find_outliers(self.gbt_s, self.taxid_2_bp)
        self.wf_data = self._find_wide_families(self.gbt_s, self.taxid_2_bp)  # List of families with lists of genera data
        self.wg_data = self._find_wide_genera(self.gbt_s, self.taxid_2_bp)

    def get_taxid_2_bp(self, gbt):
        taxid_2_bp = dict()

        def _t2b(taxon, top_taxa):
            if taxon.rank == 'genus':
                min_s = params.genus_min_sequences
                iqr_median = params.genus_iqr_median
            else:
                min_s = params.family_min_sequences
                iqr_median = params.family_iqr_median
            if len(taxon.sequence_data) >= min_s:
                bp = _boxplot_stats(taxon.sequence_data)
                bp['wide'] = (bp['iqr_median'] > iqr_median)
                bp['taxon_name'] = taxon.name  # For reports
                taxid_2_bp[taxon.taxid] = dict(taxon=taxon, top_taxa=top_taxa.name, data=bp)
                for t in taxon.children:
                    _t2b(t, top_taxa)

        for taxon in gbt.get_top_taxa():
            _t2b(taxon, taxon)
        return taxid_2_bp

    def find_outliers(self, gbt, taxid_2_bp):
        genus_taxid_2_outliers = dict()  # taxid -> list of outlier data
        for taxons in gbt.iterate_leaves():
            genus_taxon = taxons[-1]
            if not (bp := taxid_2_bp.get(genus_taxon.taxid)):
                continue
            bp_data = bp['data']
            if not (outliers := (list(bp_data['lower']) + list(bp_data['higher']))):
                continue
            #
            exp_len = int(bp_data['median'])
            res = []
            for sd in outliers:
                nucleotide = fetch_nucleotide(params, sd['taxid'])
                half = (sd['length'] + exp_len) / 2
                if sd['length'] < exp_len:
                    is_cl = lambda s: s['length'] > half
                else:
                    is_cl = lambda s: s['length'] < half
                closer_ns = [n_seq for n_seq in nucleotide if is_cl(n_seq)]
                res.append(dict(sequence=sd, exp_len=exp_len,
                                nucleotide=nucleotide, unique=len(set(s['length'] for s in nucleotide)),
                                closer_ns=closer_ns))
            genus_taxid_2_outliers[genus_taxon.taxid] = res
        return genus_taxid_2_outliers

    def _find_wide_families(self, gbt_s, taxid_2_bp):
        wf_data = []  # List of families with lists of genera data
        for taxons in gbt_s.iterate_by_rank('family'):
            f_t = taxons[-1]
            if (bp := taxid_2_bp.get(f_t.taxid)) and bp['data']['wide']:
                genera = sorted(f_t.children, key=lambda t: len(t.sequence_data), reverse=True)
                bp_d = bp['data']
                wf_data.append(dict(family=f_t.name,
                                    lengths=f_t.lengths,
                                    iqr_median=bp_d['iqr_median'],
                                    median=bp_d['median'],
                                    num_outliers=bp_d['num_lo'] + bp_d['num_hi'],
                                    genera=[dict(genus=g_t.name, lengths=g_t.lengths, taxid=g_t.taxid) for g_t in genera]))
        return sorted(wf_data, key=lambda d: d['iqr_median'], reverse=True)

    def _find_wide_genera(self, gbt_s, taxid_2_bp):
        wg = []  # Tuple (name, list of values)
        for taxons in gbt_s.iterate_leaves():
            genus_taxon = taxons[-1]
            if (bp := taxid_2_bp.get(genus_taxon.taxid)) and bp['data']['wide']:
                bp_d = bp['data']
                wg.append(dict(genus=genus_taxon.name,
                               family=taxons[-2].name,
                               lengths=genus_taxon.lengths,
                               iqr_median=bp_d['iqr_median'],
                               median=bp_d['median'],
                               num_outliers=bp_d['num_lo'] + bp_d['num_hi']))
        return sorted(wg, key=lambda d: d['iqr_median'], reverse=True)


def _boxplot_stats(sequence_data, whis=1.5):
    x = numpy.array([s['length'] for s in sequence_data])

    # medians and quartiles
    q1, med, q3 = numpy.percentile(x, [25, 50, 75])
    iqr = q3 - q1

    # get high extreme
    hival = q3 + whis * iqr
    wiskhi = x[x <= hival]
    whishi = q3 if (len(wiskhi) == 0 or numpy.max(wiskhi) < q3) else numpy.max(wiskhi)

    # get low extreme
    loval = q1 - whis * iqr
    wisklo = x[x >= loval]
    whislo = q1 if (len(wisklo) == 0 or numpy.min(wisklo) > q1) else numpy.min(wisklo)

    # compute number of outliers
    lower = [s for s in sequence_data if s['length'] < whislo]
    higher = [s for s in sequence_data if s['length'] > whishi]
    num_lo = len(lower)
    num_hi = len(higher)

    iqr_median = round(100 * iqr / med, 2)

    return dict(min_length=min(x), max_length=max(x), num_sequences=len(x),
                median=med, iqr=iqr, iqr_median=iqr_median,
                whislo=whislo, whishi=whishi,
                lower=lower, higher=higher,
                num_lo=num_lo, num_hi=num_hi,
                row=[min(x), max(x), med, iqr, whislo, whishi, num_lo, num_hi, iqr_median])


def search_biopython(db, term, retmax=18000):
    def _entrez(handle):
        data = Entrez.read(handle)
        handle.close()
        return data

    print(f'Search: {db} "{term}"')
    search = _entrez(Entrez.esearch(db=db, term=term, usehistory='y', retmax=1, rettype='json'))
    data = []
    if (count := int(search['Count'])):
        n_pages = ((count - 1) // retmax) + 1
        args = dict(db=db, query_key=search['QueryKey'], WebEnv=search['WebEnv'], retmax=retmax, rettype='json')
        for page in range(n_pages):
            print(f'    Fetch: page {page + 1} / {n_pages}, records {count}')
            data += _entrez(Entrez.esummary(**args))
            args['retstart'] = (page + 1) * retmax
    return data


def _fetch_seq_data(taxon, is_refseq, json_filename, cache_dir):
    if os.path.isfile(json_filename):
        with open(json_filename, 'r', encoding='utf-8') as _in:
            data = json.load(_in)
    else:
        q = f'txid{taxon}' if isinstance(taxon, int) or taxon.isdigit() else f'"{taxon}"'
        rs = ' AND srcdb_refseq[PROP]' if is_refseq else ''
        data = search_biopython(
            'nuccore',
            f'{q}[Organism:exp] AND "complete genome"[Title]{rs} AND (chloroplast[filter] OR plastid[filter])')

        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)
        with open(json_filename, 'w', encoding='utf-8') as _out:
            json.dump(data, _out, indent=2)

    def _species(title):
        t_fs = title.split()
        return ' '.join(t_fs[1:3] if t_fs[0].upper().startswith('UNVERIFIED') else t_fs[:2])

    return [dict(taxid=int(d['TaxId']),
                 accession=d['Caption'],
                 length=d['Length'],
                 species=_species(d['Title']),
                 title=d['Title'],
                 created=d['CreateDate'].replace('/', '-'),
                 updated=d['UpdateDate'].replace('/', '-'))
            for d in data if d['Status'] == 'live' and not d['ReplacedBy']]


def fetch_refseq(params):
    data = []
    for taxon in params.taxon:
        json_filename = os.path.join(params.cache_dir, f'{taxon}.json')
        data += _fetch_seq_data(taxon, True, json_filename, params.cache_dir)
    if params.published:
        # Note: no need for convertion in date type!
        # For valid dates, ISO formated string and date comparisons are equivalent.
        data = [d for d in data if d['updated'] <= params.published]
    return data


def fetch_nucleotide(params, taxid):
    json_filename = os.path.join(params.cache_dir, f'nucleotide_{taxid}.json')
    return _fetch_seq_data(taxid, False, json_filename, params.cache_dir)


# Box plot data
def boxplot_data(params, gbt, taxid_2_bp, taxid_2_num_genera):
    def _bd(taxon):
        n_genera = taxid_2_num_genera[taxon.taxid]
        if taxon.rank != 'genus' and n_genera == 1:  # Collaps
            return 'N', None

        n_seq = len(taxon.sequence_data)
        d = [taxon.taxid, len(taxon.children), n_genera, n_seq]
        if bp := taxid_2_bp.get(taxon.taxid):
            bp = bp['data']
            return 'C', d + bp['row'] + ['Wide' if bp['wide'] else None]
        return 'S', d + [None] * 7
    return gbt.get_taxon_data(_bd, True)


def wide_boxplot_data(params, gbt, taxid_2_bp, taxid_2_num_genera):
    def _wbd(taxon):
        if not (bp := taxid_2_bp.get(taxon.taxid)):  # No data
            return 'N', None
        bp = bp['data']
        if not bp['wide']:  # Not wide
            return 'N', None

        n_seq = len(taxon.sequence_data)
        n_genera = taxid_2_num_genera[taxon.taxid]
        return 'C', [taxon.taxid, len(taxon.children), n_genera, n_seq] + bp['row']
    return gbt.get_taxon_data(_wbd, True)


def wide_families(params, gbt, taxid_2_bp, taxid_2_num_genera):
    columns, _clade_names = gbt.column_data()
    rows = []
    last_clade_names = None

    for taxons in gbt.iterate_by_rank('family'):
        taxon = taxons[-1]
        if (bp := taxid_2_bp.get(taxon.taxid)) and bp['data']['wide']:
            clade_names, last_clade_names = remove_dupl_names(_clade_names(taxons), last_clade_names, True)
            rows.append(clade_names + [taxon.taxid, taxid_2_num_genera[taxon.taxid], len(taxon.sequence_data)] + bp['data']['row'])

            for g_taxon in sorted(taxon.children, key=(lambda t: len(t.sequence_data)), reverse=True):
                cs = _clade_names([g_taxon])
                if g_bp := taxid_2_bp.get(g_taxon.taxid):
                    g_row = g_bp['data']['row']
                else:
                    lens = g_taxon.lengths
                    if len(lens) == 1:
                        g_row = [min(lens)] + [None] * 8
                    elif len(lens) == 2:
                        g_row = [min(lens), max(lens)] + [None] * 7
                    else:
                        g_row = [min(lens), max(lens), median(lens)] + [None] * 6
                rows.append(cs + [g_taxon.taxid, None, len(g_taxon.sequence_data)] + g_row)

    return columns, rows


def wide_genera(params, gbt, taxid_2_bp):
    columns, _clade_names = gbt.column_data()
    rows = []
    last_clade_names = None
    for taxons in gbt.iterate_leaves():
        genus_taxon = taxons[-1]
        if (bp := taxid_2_bp.get(genus_taxon.taxid)) and bp['data']['wide']:
            clade_names, last_clade_names = remove_dupl_names(_clade_names(taxons), last_clade_names, True)
            rows.append(clade_names + [genus_taxon.taxid, len(genus_taxon.sequence_data)] + bp['data']['row'])
    return columns, rows


# Outliers
def genera_outliers_all(params, gbt, taxid_2_bp, genus_taxid_2_outliers):
    columns, _clade_names = gbt.column_data()
    rows = []
    last_clade_names = None
    for taxons in gbt.iterate_leaves():
        genus_taxon = taxons[-1]
        if not (outlier_data := genus_taxid_2_outliers.get(genus_taxon.taxid)):
            continue
        #
        clade_names, last_clade_names = remove_dupl_names(_clade_names(taxons), last_clade_names, True)
        bp = taxid_2_bp[genus_taxon.taxid]['data']
        wide = 'Wide' if bp.get('wide') else None
        for od in outlier_data:
            sd = od['sequence']
            rows.append(clade_names +
                        [sd[a] for a in ('species', 'taxid', 'accession', 'updated', 'length')] +
                        [od['exp_len'], bp['iqr_median'], wide, len(od['nucleotide']), od['unique'], 'YES' if od['closer_ns'] else None])
            clade_names = [None] * len(clade_names)
    return columns, rows


def genera_outliers(params, gbt, taxid_2_bp, genus_taxid_2_outliers):
    columns, _clade_names = gbt.column_data()
    rows = []
    last_clade_names = None
    for taxons in gbt.iterate_leaves():
        genus_taxon = taxons[-1]
        if not (outlier_data := genus_taxid_2_outliers.get(genus_taxon.taxid)):
            continue
        if not any(od['closer_ns'] for od in outlier_data):
            continue
        #
        clade_names, last_clade_names = remove_dupl_names(_clade_names(taxons), last_clade_names, True)
        bp = taxid_2_bp[genus_taxon.taxid]['data']
        wide = 'Wide' if bp.get('wide') else None
        for od in outlier_data:
            sd = od['sequence']
            seq_row = clade_names + \
                [sd[a] for a in ('species', 'taxid', 'accession', 'updated', 'length')] + \
                [od['exp_len'], bp['iqr_median'], wide, len(od['nucleotide']), od['unique'], len(od['closer_ns'])]
            for cl_ns in od['closer_ns']:
                rows.append(seq_row + [cl_ns['accession'], cl_ns['length'], cl_ns['updated']])
                seq_row = [None] * len(seq_row)
    return columns, rows


#
def export_excel(filename, params, results):
    sheets = []

    # Sheet RefSeq data
    cs = ['species', 'taxid', 'accession', 'length', 'created', 'updated', 'title']
    columns, rows = results.gbt_l.get_table_data((lambda sd: [sd[a] for a in cs]), (lambda sd: sd['species']), True)
    sheets.append(('RefSeq', columns + [c.capitalize() for c in cs], rows))

    # Boxplot data
    columns, rows = boxplot_data(params, results.gbt_s, results.taxid_2_bp, results.taxid_2_num_genera)
    cs_nums = ['Taxid', 'Num Lower', 'Num Genera', 'Num Sequences']
    cq_bp = ['Min Length', 'Max Length', 'Median', 'IQR', 'WhisLo', 'WhisHi', 'NumLo', 'NumHi', 'IQR/Median']
    sheets.append(('Distributions', columns + cs_nums + cq_bp + ['Wide'], rows))

    columns, rows = wide_boxplot_data(params, results.gbt_s, results.taxid_2_bp, results.taxid_2_num_genera)
    sheets.append(('WideAll', columns + cs_nums + cq_bp, rows))

    columns, rows = wide_families(params, results.gbt_s, results.taxid_2_bp, results.taxid_2_num_genera)
    sheets.append(('WideFamilies', columns + ['Taxid', 'Num Genera', 'Num Sequences'] + cq_bp, rows))

    columns, rows = wide_genera(params, results.gbt_s, results.taxid_2_bp)
    sheets.append(('WideGenera', columns + ['Taxid', 'Num Sequences'] + cq_bp, rows))

    # Outliers
    columns, rows = genera_outliers_all(params, results.gbt_s, results.taxid_2_bp, results.genus_taxid_2_outliers)
    cs = ['Species', 'Taxid', 'Accession', 'Published', 'Length', 'ExpLength', 'IQR/Median', 'WideGenus', 'Nucleotide', 'Unique', 'Closer']
    sheets.append(('Outliers', columns + cs, rows))

    columns, rows = genera_outliers(params, results.gbt_s, results.taxid_2_bp, results.genus_taxid_2_outliers)
    cs += ['CloserSeq', 'CloserLength', 'CloserPublished']
    sheets.append(('Alternatives', columns + cs, rows))

    sheets_to_excel(filename, sheets)


#
def create_summary(filename, params, results):
    gbt = results.gbt_s
    nsb = _NumSpeciesBelow(params.cache_dir)
    clade_2_num_seq = dict((t.name, len(t.sequence_data)) for t in gbt.get_top_taxa())
    _perc = lambda v: f'{round(100 * v, 2)}%'

    def _pp(num_all, num_part):
        return f'{num_part} ({_perc(num_part / num_all)})'

    def _bp(rank, label, min_seq, iqr_median):
        bps = [bp['data'] for bp in results.taxid_2_bp.values() if bp['taxon'].rank == rank]

        in_taxa = sum(bp['num_sequences'] for bp in bps)
        ncbi_in_taxa = 0
        if params.clade:
            clade_bps = dict((t.name, []) for t in gbt.get_top_taxa())  # Use tree sorting
            clade_min_species = dict((t.name, (10000000, ''))for t in gbt.get_top_taxa())
            clade_max_species = dict((t.name, (0, ''))for t in gbt.get_top_taxa())

        for bp in results.taxid_2_bp.values():
            if bp['taxon'].rank == rank:
                ns = (nsb.get(bp['taxon'].taxid), bp['taxon'].name)
                ncbi_in_taxa += ns[0]
                if params.clade:
                    clade = bp['top_taxa']
                    clade_bps[clade].append(bp['data'])
                    clade_min_species[clade] = min(clade_min_species[clade], ns)
                    clade_max_species[clade] = max(clade_max_species[clade], ns)

        #
        iqrm_ps = [0.3, 0.5] + list(range(1, 11))
        iqrm_in_clases = ''

        #
        if params.clade:
            clade_bps = list(clade_bps.items())  # Shorter
            taxa_in_clades = ' ; ' + ', '.join(f'{c} - {len(_bps)}' for c, _bps in clade_bps)
            seq_in_clades = ' ; ' + ', '.join(f'{c} - {sum(bp["num_sequences"] for bp in _bps)}' for c, _bps in clade_bps)
            min_in_clades = [(c, min(_bps, key=lambda x: x['num_sequences'])) for c, _bps in clade_bps]
            min_in_clades = ', '.join(f'{c} - {bp["taxon_name"]} ({bp["num_sequences"]})' for c, bp in min_in_clades)
            max_in_clades = [(c, max(_bps, key=lambda x: x['num_sequences'])) for c, _bps in clade_bps]
            max_in_clades = ', '.join(f'{c} - {bp["taxon_name"]} ({bp["num_sequences"]})' for c, bp in max_in_clades)
            wide_in_clades = ' ; ' + ', '.join(f'{c} - {sum(int(bp["wide"]) for bp in _bps)}' for c, _bps in clade_bps)
            taxa_outliers_in_clades = ' ; ' + ', '.join(f'{c} - {_pp(len(_bps), sum(int(bool(bp["num_lo"] or bp["num_hi"])) for bp in _bps))}' for c, _bps in clade_bps)
            outliers_in_clades = ' ; ' + ', '.join(f'{c} - {_pp(clade_2_num_seq[c], sum(bp["num_lo"] + bp["num_hi"] for bp in _bps))}' for c, _bps in clade_bps)
            for c, _bps in clade_bps:
                _ims = [x['iqr_median'] for x in _bps]
                iqrm_in_clases += f"\n    {c:>21}  : {', '.join(f'{p}-{sum(int(x > p) for x in _ims)}' for p in iqrm_ps)}"
            wide_genus_in_clades = ' ; ' + ', '.join(f'{c} - {_pp(len(_bps), sum(int(bp["iqr_median"] > params.genus_iqr_median) for bp in _bps))}' for c, _bps in clade_bps)
            wide_family_in_clades = ' ; ' + ', '.join(f'{c} - {_pp(len(_bps), sum(int(bp["iqr_median"] > params.family_iqr_median) for bp in _bps))}' for c, _bps in clade_bps)

            min_max_clade = f"""
    minimal sequences      : {min_in_clades}
    maximal sequences      : {max_in_clades}"""
            min_max_ncbi_clade = f"""
    minimal NCBI species   : {', '.join(f'{c} - {t} ({n})' for c, (n, t) in clade_min_species.items())}
    maximal NCBI species   : {', '.join(f'{c} - {t} ({n})' for c, (n, t) in clade_max_species.items())}"""

        else:
            taxa_in_clades = ''
            seq_in_clades = ''
            wide_in_clades = ''
            taxa_outliers_in_clades = ''
            outliers_in_clades = ''
            wide_genus_in_clades = ''
            wide_family_in_clades = ''
            min_max_clade = ''
            min_max_ncbi_clade = ''

        iqr_meds = [x['iqr_median'] for x in bps]

        num_los = sum(bp['num_lo'] for bp in bps)
        num_his = sum(bp['num_hi'] for bp in bps)
        num_all = num_los + num_his
        num_taxa = sum(int(bool(bp['num_lo'] or bp['num_hi'])) for bp in bps)

        return f"""{label}
  Data
    with enough sequences  : {len(bps):>4}  (minimal sequences {min_seq}){taxa_in_clades}
    sequences in all calcs : {in_taxa:>4}  ({_perc(in_taxa / len(results.sequence_data))}){seq_in_clades}{min_max_clade}
    NCBI species in taxa   : {ncbi_in_taxa}  ({_perc(in_taxa / ncbi_in_taxa)}){min_max_ncbi_clade}
  Box plot
    wide distribution      : {sum(int(x > iqr_median) for x in iqr_meds):>4}{wide_in_clades}
         by genus perc     : {_pp(len(bps), sum(int(x > params.genus_iqr_median) for x in iqr_meds))}{wide_genus_in_clades}
         by family  perc   : {_pp(len(bps), sum(int(x > params.family_iqr_median) for x in iqr_meds))}{wide_family_in_clades}
    num IRQ/media > 1-10%  : {', '.join(f'{p}-{sum(int(x > p) for x in iqr_meds)}' for p in iqrm_ps)}{iqrm_in_clases}
    num outlier taxa       : {_pp(len(bps), num_taxa)}{taxa_outliers_in_clades}
    num outlier seqs.      : {_pp(in_taxa, num_all):>4}  (lower {num_los}, higher {num_his}){outliers_in_clades}
"""

    #
    def _g_out(label, test):
        gs = [g for g_t, g in results.genus_taxid_2_outliers.items() if test(g_t)]
        return f"""{label}
  num genera               : {len(gs):>5}
  species to test          : {sum(len(g) for g in gs):>5}
  nucleotide downloaded    : {sum(sum(len(s['nucleotide']) for s in g) for g in gs):>5}
  species with alternative : {sum(sum(int(s['unique'] > 1) for s in g) for g in gs):>5}
  with better alternative  : {sum(sum(int(bool(s['closer_ns'])) for s in g) for g in gs):>5}
  newer sequences          : {sum(sum(int(any(s['sequence']['updated'] < c['updated'] for c in s['closer_ns'])) for s in g) for g in gs):>5}
"""

    #
    seq_in_clades = ' ; ' + ', '.join(f'{t.name} - {len(t.sequence_data):>4}' for t in gbt.get_top_taxa())
    gen_in_clades = ' ; ' + ', '.join(f'{t.name} - {t.num_of_rank("genus"):>4}' for t in gbt.get_top_taxa())
    fam_in_clades = ' ; ' + ', '.join(f'{t.name} - {t.num_of_rank("family"):>4}' for t in gbt.get_top_taxa())

    wfs = [f'{"Family":<15}{"Genera":>10}{"Seq":>10}{"IQR/median":>10}{"Outliers":>10}{"MoreDist":>10}   LowerWide']
    wf_by_num = defaultdict(int)
    for wf in results.wf_data:
        min_s = min(params.genus_min_sequences, len(wf['lengths']) // 10)
        g_ms = [median(g['lengths']) for g in wf['genera'] if len(g['lengths']) >= min_s]
        more_dist = g_ms and (100 * (max(g_ms) - min(g_ms)) / wf['median']) > params.genus_iqr_median
        lower = [g['genus'] for g in wf['genera'] if (bp := results.taxid_2_bp.get(g['taxid'])) and bp['data']['wide']]
        wfs.append(f'{wf["family"]:<15}{len(wf["genera"]):>10}{len(wf["lengths"]):>10}{wf["iqr_median"]:>10}{"+" if wf["num_outliers"] else "-":>10}{"+" if more_dist else "-":>10}   {", ".join(lower)}')
        wf_by_num[int(more_dist) + int(bool(wf["num_outliers"])) + int(bool(lower))] += 1
    wfs = '  Described\n' + '\n'.join(('    ' + line) for line in wfs)
    wfs += '\n  By num properties: ' + '; '.join(f'{n} - {c}' for n, c in sorted(wf_by_num.items())) + '\n\n'

    wgs = [f'{"Genus":<15}{"Family":<15}{"IQR/median":>10}{"Seq":>10}{"Outliers":>10}']
    for wg in results.wg_data:
        wgs.append(f'{wg["genus"]:<15}{wg["family"]:<15}{wg["iqr_median"]:>10}{len(wg["lengths"]):>10}{"+" if wg["num_outliers"] else "-":>10}')
    wgs = '  Described\n' + '\n'.join(('    ' + line) for line in wgs)

    #
    outliers_1 = []
    for outliers in results.genus_taxid_2_outliers.values():
        for outlier in outliers:
            if outlier['closer_ns']:
                # s_l = outlier['sequence']['length']
                exp_len = outlier['exp_len']
                s_exp = abs(outlier['sequence']['length'] - exp_len)
                closests = min(outlier['closer_ns'], key=lambda c: abs(c['length'] - exp_len) / s_exp)
                outliers_1.append((abs(closests['length'] - exp_len) / s_exp, exp_len, outlier['sequence'], closests))
    out_top = [f'  {"Species":<20}{"Accession":>12}{"Published":>12}{"Length":>10}{"Median":>10}{"ClLength":>10}{"Closer":>10}{"ClPublished":>12}']
    for _, gen_med, seq, cl in sorted(outliers_1)[:10]:
        out_top.append(f'  {seq["species"]:<20}{seq["accession"]:>12}{seq["updated"]:>12}{seq["length"]:>10}{gen_med:>10}{cl["length"]:>10}{cl["accession"]:>10}{cl["updated"]:>12}')
    out_top = '\n'.join(out_top)

    #
    summary = f"""
GENERAL

Number of sequences    : {len(results.sequence_data):>5}{seq_in_clades}
Number of genera       : {gbt.num_of_rank('genus'):>5}{gen_in_clades}
Number of families     : {gbt.num_of_rank('family'):>5}{fam_in_clades}


BOX PLOT DATA

{_bp('family', 'Families', params.family_min_sequences, params.family_iqr_median)}{wfs}
{_bp('genus', 'Genera', params.genus_min_sequences, params.genus_iqr_median)}{wgs}

GENERA OUTLIERS

{_g_out('All', lambda gt: True)}
{_g_out('Wide', lambda gt: results.taxid_2_bp[gt]['data'].get('wide'))}
{_g_out('Narrow', lambda gt: not results.taxid_2_bp[gt]['data'].get('wide'))}

Alterantive sample
{out_top}

"""

    #
    print(summary)
    with open(filename, 'w') as _out:
        _out.write(summary)
    nsb.save()


#
def create_figure_data(base_name, fd_filename, params, results):
    # Data, list of dicts with family data:
    #   clade, order, family, num_species, num_sequences, list of publishing dates, list of lengths
    nsb = _NumSpeciesBelow(params.cache_dir)
    if params.clade:
        c_m = lambda taxons: taxons[0].name
        o_idx, f_idx = 1, 2
    else:
        c_m = lambda taxons: ''
        o_idx, f_idx = 0, 1
    data = [dict(clade=c_m(taxons),
                 order=taxons[o_idx].name,
                 family=taxons[f_idx].name,
                 num_species=nsb.get(taxons[f_idx].taxid),
                 num_sequences=len(taxons[f_idx].sequence_data),
                 dates=[s['updated'] for s in taxons[f_idx].sequence_data],
                 lengths=taxons[f_idx].lengths)
            for taxons in results.gbt_s.iterate_by_rank('family')
            if len(taxons[-1].sequence_data) >= params.family_min_sequences]
    nsb.save()

    # IQR/median values
    iqr_median = dict(families=[bp['data']['iqr_median'] for bp in results.taxid_2_bp.values() if bp['taxon'].rank == 'family'],
                      genera=[bp['data']['iqr_median'] for bp in results.taxid_2_bp.values() if bp['taxon'].rank == 'genus'])

    # Clades
    clade_names = dict((top_t.name, [f_t.name for f_t in f_ts]) for top_t, f_ts in results.gbt_s.group_by_rank('family'))
    for top_t, g_ts in results.gbt_s.group_by_rank('genus'):
        clade_names[top_t.name] += [g_t.name for g_t in g_ts]

    #
    output_data = dict(base_name=base_name, data=data, iqr_median=iqr_median,
                       wf=results.wf_data, wg=results.wg_data, clade_names=clade_names)
    with open(fd_filename, 'w', encoding='utf-8') as _out:
        json.dump(output_data, _out, indent=2)

#
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Survey NCBI cpDNA data for relation between species closeness and sequence lengths')

    parser.add_argument('taxon', nargs='+', help='Taxa name to collect data')

    parser.add_argument('-g', '--genus-min-sequences', type=int, default=10,
                        help='Minimal number of sequences to use for genus statistics')
    parser.add_argument('-f', '--family-min-sequences', type=int, default=20,
                        help='Minimal number of sequences to use for family statistics')
    parser.add_argument('-G', '--genus-iqr-median', type=float, default=1.,
                        help='IQR/median ratio (in %) to take genus distribution as wide')
    parser.add_argument('-F', '--family-iqr-median', type=float, default=10.,
                        help='IQR/median ratio (in %) to take genus distribution as wide')

    # Additional arguments
    parser.add_argument('-p', '--published', help='Max date of update. In ISO format yyyy-mm-dd.')
    parser.add_argument('-c', '--clade', action='append', help='Clade(s) to use for grouping')
    parser.add_argument('-E', '--email', help='Email address, for Entrez')
    parser.add_argument('-C', '--cache-dir', default='cache_data', help='Cache directory for NCBI queries')

    parser.add_argument('-i', '--images', action='store_true', help='Create images')
    parser.add_argument('--dpi', type=int, help='DPI value use for saving images')
    parser.add_argument('--svg', action='store_true', help='Store image in SVG format')

    params = parser.parse_args()

    #
    if params.email:
        Entrez.email = params.email

    results = _Results(params)
    base_name = f'{"_".join(sorted(params.taxon))}_g{params.genus_min_sequences}_f{params.family_min_sequences}'
    export_excel(f'{base_name}_results.xlsx', params, results)
    create_summary(f'{base_name}_summary.txt', params, results)
    fd_filename = f'{base_name}_figures_data.json'
    create_figure_data(base_name, fd_filename, params, results)
    #
    if params.images:
        exe = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'survey_figures.py')
        cmd_args = ['--no-preview']
        if params.dpi:
            cmd_args.extend(['--dpi', str(params.dpi)])
        if params.svg:
            cmd_args.append('--svg')
        for im_type in ('data', 'im', 'wf', 'wg'):
            subprocess.run([exe, fd_filename, im_type] + cmd_args, check=True)
