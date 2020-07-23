'''The API calls'''

from io import BytesIO
import json
from flask import (
    abort,
    g,
    jsonify,
    redirect,
    request,
    Response,
    send_file,
    stream_with_context,
)
import re
import sqlalchemy
from sqlalchemy import (
    cast,
    desc as sql_desc,
    distinct,
    Float,
    func,
)
from . import app, taxtree
from .search import (
    core_search,
    format_results,
    json_stats,
    available_term_by_category,
)
from .search.helpers import UnknownQueryError
from .search_parser import Query
from .models import (
    db,
    BgcType,
    Region,
    DnaSequence,
    Filename,
    Genome,
    Taxa,
    t_rel_regions_types,
)
from .errors import TooManyResults
from .legacy import dbv1_accessions


MIME_TYPE_MAP = {
    'json': 'application/json',
    'csv': 'text/csv',
    'fasta': 'application/fasta',
    'fastaa': 'application/fasta',
}

FASTA_LIMITS = {
    'cluster': 100,
    'gene': 500,
    'domain': 500
}

SAFE_IDENTIFIER_PATTERN = re.compile('[^A-Za-z0-9_.]+', re.UNICODE)


@app.route('/api/v1.0/version')
def get_version():
    '''display the API version'''
    from .version import __version__ as api_version
    ret = {
        'api': api_version,
        'sqlalchemy': sqlalchemy.__version__
    }
    return jsonify(ret)


def _common_stats():
    """Get the stats shared by the v1 and v2 version of the call"""

    num_clusters = Region.query.filter(Region.minimal.is_(False)).count()

    num_genomes = Genome.query.count()

    num_sequences = DnaSequence.query.count()

    clusters = []
    sub = db.session.query(t_rel_regions_types.c.bgc_type_id, func.count(1).label('count')) \
                    .join(Region).filter(Region.minimal.is_(False)) \
                    .group_by(t_rel_regions_types.c.bgc_type_id).subquery()
    ret = db.session.query(BgcType.term, BgcType.description, sub.c.count).join(sub) \
                    .order_by(sub.c.count.desc(), BgcType.term)
    for cluster in ret:
        clusters.append({'name': cluster.term, 'description': cluster.description, 'count': cluster.count})

    ret = db.session.query(Taxa.tax_id, Taxa.genus, Taxa.species, func.count(DnaSequence.accession).label('tax_count')) \
                    .join(Genome, Genome.tax_id == Taxa.tax_id).join(DnaSequence, DnaSequence.genome_id == Genome.genome_id) \
                    .group_by(Taxa.tax_id).order_by(sql_desc('tax_count')).limit(1).first()
    top_seq_taxon = ret.tax_id
    top_seq_species = '{r.genus} {r.species}'.format(r=ret)
    top_seq_taxon_count = ret.tax_count

    stats = {
        'num_clusters': num_clusters,
        'num_genomes': num_genomes,
        'num_sequences': num_sequences,
        'top_seq_taxon': top_seq_taxon,
        'top_seq_taxon_count': top_seq_taxon_count,
        'top_seq_species': top_seq_species,
        'clusters': clusters,
    }
    return stats


@app.route('/api/v1.0/stats')
def get_stats_v1():
    '''contents for the stats page'''
    stats = _common_stats()

    ret = db.session.query(Taxa.tax_id, Taxa.genus, Taxa.species, Taxa.strain,
                           DnaSequence.acc,
                           func.count(distinct(Region.region_number)).label('bgc_count'),
                           func.count(distinct(DnaSequence.acc)).label('seq_count'),
                           (cast(func.count(distinct(Region.region_number)), Float) / func.count(distinct(DnaSequence.acc))).label('clusters_per_seq')) \
                    .join(Genome).join(DnaSequence).join(Region) \
                    .group_by(Taxa.tax_id, DnaSequence.acc).order_by(sql_desc('clusters_per_seq')).limit(1).first()
    stats['top_secmet_taxon'] = ret.tax_id
    stats['top_secmet_species'] = '{r.genus} {r.species} {r.strain}'.format(r=ret)
    stats['top_secmet_acc'] = ret.acc
    stats['top_secmet_taxon_count'] = ret.clusters_per_seq

    return jsonify(stats)


@app.route('/api/v2.0/stats')
def get_stats_v2():
    """contents for the stats page"""
    stats = _common_stats()

    ret = db.session.query(Taxa.tax_id, Taxa.genus, Taxa.species, Taxa.strain,
                           Genome.assembly_id,
                           func.count(distinct(Region.region_number)).label('bgc_count'),
                           func.count(distinct(Genome.assembly_id)).label('seq_count'),
                           (cast(func.count(distinct(Region.region_number)), Float) / func.count(distinct(Genome.assembly_id))).label('clusters_per_seq')) \
                    .join(Genome).join(DnaSequence).join(Region) \
                    .filter(Genome.assembly_id is not None).filter(Region.minimal.is_(False)) \
                    .group_by(Taxa.tax_id, Genome.assembly_id).order_by(sql_desc('clusters_per_seq')).limit(1).first()
    stats['top_secmet_taxon'] = ret.tax_id
    stats['top_secmet_species'] = '{r.genus} {r.species} {r.strain}'.format(r=ret)
    stats['top_secmet_assembly_id'] = ret.assembly_id
    stats['top_secmet_taxon_count'] = ret.bgc_count

    return jsonify(stats)


@app.route('/api/v1.0/tree/secmet')
def get_sec_met_tree():
    '''Get the jsTree structure for secondary metabolite clusters'''
    ret = db.session.query(Region.region_id, Region.region_number,
                           DnaSequence.accession,
                           BgcType.term, BgcType.description,
                           Taxa.genus, Taxa.species, Taxa.strain, Genome.assembly_id) \
                    .join(DnaSequence, DnaSequence.accession == Region.record_accession).join(Genome).join(Taxa) \
                    .join(t_rel_regions_types).join(BgcType) \
                    .order_by(BgcType.description, Taxa.genus, Taxa.species, DnaSequence.accession, Region.region_number) \
                    .all()

    clusters = []
    types = {}

    for entry in ret:
        types[entry.term] = entry.description
        species = entry.species if entry.species != 'Unclassified' else 'sp.'
        assembly_id = entry.assembly_id.split('.')[0] if entry.assembly_id else None
        name = '{} {} {}'.format(entry.genus, species, entry.strain)
        clusters.append({
            "id": "{}_c{}_{}".format(entry.accession, entry.region_number, entry.term),
            "parent": entry.term,
            "text": "{} {} Region {}".format(name, entry.accession, entry.region_number),
            "assembly_id": assembly_id,
            "region_number": entry.region_number,
            "type": "cluster",
        })

    for name, desc in sorted(list(types.items()), reverse=True):
        clusters.insert(0, {
            "id": name,
            "parent": "#",
            "text": desc,
            "state": {
                "disabled": True
            }
        })

    tree = clusters
    return jsonify(tree)


def _get_taxon_tree_node(tree_id):
    HANDLERS = {
        'superkingdom': taxtree.get_phylum,
        'phylum': taxtree.get_class,
        'class': taxtree.get_order,
        'order': taxtree.get_family,
        'family': taxtree.get_genus,
        'genus': taxtree.get_species,
        'species': taxtree.get_strains,
    }

    if tree_id == '1':
        tree = taxtree.get_superkingdom()
    else:
        params = tree_id.split('_')
        taxlevel = params[0]
        params = params[1:]
        handler = HANDLERS.get(taxlevel, lambda x: [])
        tree = handler(params)

    return tree


@app.route('/api/v1.0/tree/taxa')
def get_taxon_tree():
    '''Get the jsTree structure for all taxa'''
    tree_id = request.args.get('id', '1')
    tree = _get_taxon_tree_node(tree_id)

    return jsonify(tree)


@app.route('/api/v1.0/tree/taxa/massload')
def get_taxon_tree_massload():
    tree_ids = request.args.get('id', '1')
    id_list = tree_ids.split(',')

    multitree = {}
    for tree_id in id_list:
        multitree[tree_id] = _get_taxon_tree_node(tree_id)

    return jsonify(multitree)


@app.route('/api/v1.0/tree/taxa/search')
def search_taxon_tree():
    search = request.args.get('str', None)
    if not search:
        return jsonify([])

    search_path = taxtree.search(search)

    return jsonify(search_path)


@app.route('/api/v1.0/search', methods=['POST'])
def search():
    try:
        if 'query' not in request.json:
            query = Query.from_string(request.json.get('search_string', ''))
        else:
            query = Query.from_json(request.json['query'])
    except ValueError:
        abort(400)

    if query.return_type != 'json':
        abort(400)

    try:
        offset = int(request.json.get('offset', '0'))
    except ValueError:
        offset = 0

    try:
        paginate = int(request.json.get('paginate', '50'))
    except ValueError:
        paginate = 50

    try:
        clusters = format_results(query, core_search(query))
    except UnknownQueryError:
        abort(400)
    stats = json_stats(clusters)

    total = len(clusters)

    if paginate > 0:
        end = min(offset + paginate, total)
    else:
        end = total

    result = {
        'total': total,
        'clusters': clusters[offset:end],
        'offset': offset,
        'paginate': paginate,
        'stats': stats,
    }

    return jsonify(result)


@app.route('/api/v1.0/export', methods=['POST'])
def export():
    '''Export the search results as CSV file'''
    try:
        if 'query' not in request.json:
            query = Query.from_string(request.json.get('search_string', ''), return_type='csv')
        else:
            query = Query.from_json(request.json['query'])
    except ValueError:
        abort(400)

    try:
        offset = int(request.json.get('offset', '0'))
    except ValueError:
        offset = 0

    try:
        paginate = int(request.json.get('paginate', '0'))
    except ValueError:
        paginate = 0

    return_type = query.return_type
    search_type = query.search_type

    if return_type not in ('json', 'csv', 'fasta', 'fastaa'):
        abort(400)

    try:
        search_results = core_search(query)
    except UnknownQueryError:
        abort(400)
    g.verbose = query.verbose
    if query.verbose:
        g.search_str = str(query)

    total = len(search_results)

    if paginate > 0:
        end = min(offset + paginate, total)
    else:
        end = total

    search_results = search_results[offset:end]

    limit = FASTA_LIMITS.get(search_type, 100)

    if return_type.startswith('fasta') and len(search_results) > limit:
        raise TooManyResults('More than {limit} search results for FASTA {search} download ({number} found), please specify a smaller query.'.format(
            limit=limit, search=search_type, number=len(search_results)))

    found_bgcs = format_results(query, search_results)
    filename = 'asdb_search_results.{}'.format(return_type)
    if query.return_type == 'json':
        found_bgcs = [json.dumps(found_bgcs)]

    handle = BytesIO()
    for line in found_bgcs:
        handle.write('{}\n'.format(line).encode('utf-8'))

    handle.seek(0)

    mime_type = MIME_TYPE_MAP.get(query.return_type, None)

    return send_file(handle, mimetype=mime_type, attachment_filename=filename, as_attachment=True)


@app.route('/api/v1.0/export/<search_type>/<return_type>')
def export_get(search_type, return_type):
    '''Export the search results as a file'''

    search_string = request.args.get('search', '')
    if search_string == '':
        abort(400)

    if return_type not in ('json', 'csv', 'fasta', 'fastaa'):
        abort(400)

    query = Query.from_string(search_string, search_type=search_type, return_type=return_type)

    g.verbose = False
    try:
        search_results = core_search(query)
    except UnknownQueryError:
        abort(400)
    if len(search_results) > 100 and search_type == 'cluster' and return_type == 'fasta':
        raise TooManyResults('More than 100 search results for FASTA cluster download, please specify a smaller query.')
    found_bgcs = format_results(query, search_results)
    if query.return_type == 'json':
        found_bgcs = [json.dumps(found_bgcs)]

    def generate():
        for line in found_bgcs:
            yield line + '\n'

    mime_type = MIME_TYPE_MAP.get(query.return_type, None)

    return Response(stream_with_context(generate()), mimetype=mime_type)


@app.route('/api/v1.0/genome/<identifier>')
def show_genome(identifier):
    '''show information for a genome by identifier'''
    query = Query.from_string('[acc]{}'.format(identifier))
    found_bgcs = format_results(query, core_search(query))

    return jsonify(found_bgcs)


@app.route('/api/v1.0/assembly/<identifier>')
def show_assembly(identifier):
    """show information for an assembly by identifier"""
    query = Query.from_string('[assembly]{}'.format(identifier))
    found_bgcs = format_results(query, core_search(query))

    return jsonify(found_bgcs)


@app.route('/api/v1.0/available/<category>/<term>')
def list_available(category, term):
    '''list available terms for a given category'''
    return jsonify(available_term_by_category(category, term))


def _canonical_assembly_id(identifier):
    """Turn the identifier into an ID usable for lookups."""
    safe_id = SAFE_IDENTIFIER_PATTERN.sub('', identifier).split('.')[0]
    if safe_id in dbv1_accessions:
        return safe_id, True

    res = db.session.query(Genome.assembly_id).filter(Genome.assembly_id.ilike("{}%".format(safe_id))).first()
    if res:
        return res.assembly_id, False

    res = db.session.query(Genome.assembly_id) \
                    .join(DnaSequence) \
                    .filter(DnaSequence.acc.ilike("{}%".format(safe_id))) \
                    .first()
    if res:
        return res.assembly_id, False

    abort(404)


@app.route('/api/v1.0/goto/<identifier>')
@app.route('/go/<identifier>')
def goto(identifier):
    safe_id, is_v1 = _canonical_assembly_id(identifier)
    if is_v1:
        return redirect("https://antismash-dbv1.secondarymetabolites.org/output/{}/index.html".format(safe_id))

    return redirect("/output/{}/index.html".format(safe_id))


@app.route('/api/v1.0/goto/<identifier>/cluster/<int:number>')
@app.route('/go/<identifier>/<int:number>')
def goto_cluster(identifier, number):
    safe_id, is_v1 = _canonical_assembly_id(identifier)
    if is_v1:
        return redirect("https://antismash-dbv1.secondarymetabolites.org/output/{}/index.html#cluster-{}".format(safe_id, number))

    return redirect("/output/{}/index.html#cluster-{}".format(safe_id, number))


def _get_base_url(identifier):
    safe_id = SAFE_IDENTIFIER_PATTERN.sub('', identifier).split('.')[0]
    ret = db.session.query(Filename.assembly_id, Filename.base_filename) \
            .filter(Filename.assembly_id == safe_id).first()
    if not ret:
        abort(404)

    return "/output/{r.assembly_id}/{r.base_filename}".format(r=ret)


@app.route('/api/v1.0/download/genbank/<identifier>')
def download_genbank(identifier):
    url = _get_base_url(identifier)
    return redirect("{}.final.gbk".format(url))


@app.route('/api/v1.0/download/table/<identifier>')
def download_table(identifier):
    url = _get_base_url(identifier)
    return redirect("{}.geneclusters.xls".format(url))


@app.route('/api/v1.0/download/genbank/<identifier>/cluster/<int:number>')
def download_cluster(identifier, number):
    url = _get_base_url(identifier)
    return redirect("{}.cluster{:03d}.gbk".format(url, number))
