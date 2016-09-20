'''The API calls'''

import StringIO
from flask import (
    jsonify,
    request,
    send_file,
)
import sqlalchemy
from sqlalchemy import (
    cast,
    desc as sql_desc,
    distinct,
    Float,
    func,
)
from . import app
from .helpers import get_db
from .search import (
    search_bgcs,
    create_cluster_csv,
    available_term_by_category,
)
from .models import (
    db,
    BgcType,
    BiosyntheticGeneCluster as Bgc,
    DnaSequence,
    Genome,
    Locus,
    Taxa,
    t_rel_clusters_types,
)
import sql


@app.route('/api/v1.0/version')
def get_version():
    '''display the API version'''
    from .version import __version__ as api_version
    ret = {
        'api': api_version,
        'sqlalchemy': sqlalchemy.__version__
    }
    return jsonify(ret)


@app.route('/api/v1.0/stats')
def get_stats():
    '''contents for the stats page'''
    num_clusters = Bgc.query.count()

    num_genomes = Genome.query.count()

    num_sequences = DnaSequence.query.count()

    clusters = []

    sub = db.session.query(t_rel_clusters_types.c.bgc_type_id, func.count(1).label('count')) \
                    .group_by(t_rel_clusters_types.c.bgc_type_id).subquery()
    ret = db.session.query(BgcType.term, BgcType.description, sub.c.count).join(sub) \
                    .order_by(sub.c.count.desc())
    for cluster in ret:
        clusters.append({'name': cluster.term, 'description': cluster.description, 'count': cluster.count})

    ret = db.session.query(Taxa.tax_id, Taxa.genus, Taxa.species, func.count(DnaSequence.acc).label('tax_count')) \
                    .join(Genome).join(DnaSequence) \
                    .group_by(Taxa.tax_id).order_by(sql_desc('tax_count')).limit(1).first()
    top_seq_taxon = ret.tax_id
    top_seq_taxon_count = ret.tax_count


    ret = db.session.query(Taxa.tax_id, Taxa.genus, Taxa.species, func.count(distinct(Bgc.bgc_id)).label('bgc_count'),
                           func.count(distinct(DnaSequence.acc)).label('seq_count'),
                           (cast(func.count(distinct(Bgc.bgc_id)), Float) / func.count(distinct(DnaSequence.acc))).label('clusters_per_seq')) \
                    .join(Genome).join(DnaSequence).join(Locus).join(Bgc) \
                    .group_by(Taxa.tax_id).order_by(sql_desc('clusters_per_seq')).limit(1).first()
    top_secmet_taxon = ret.tax_id
    top_secmet_species = ret.species
    top_secmet_taxon_count = ret.clusters_per_seq

    stats = {
        'num_clusters': num_clusters,
        'num_genomes': num_genomes,
        'num_sequences': num_sequences,
        'top_seq_taxon': top_seq_taxon,
        'top_seq_taxon_count': top_seq_taxon_count,
        'top_secmet_taxon': top_secmet_taxon,
        'top_secmet_taxon_count': top_secmet_taxon_count,
        'top_secmet_species': top_secmet_species,
        'clusters': clusters,
    }

    return jsonify(stats)


@app.route('/api/v1.0/tree/secmet')
def get_sec_met_tree():
    '''Get the jsTree structure for secondary metabolite clusters'''
    ret = db.session.query(Bgc.bgc_id, Bgc.cluster_number,
                           DnaSequence.acc,
                           BgcType.term, BgcType.description,
                           Taxa.genus, Taxa.species) \
                    .join(t_rel_clusters_types).join(BgcType).join(Locus) \
                    .join(DnaSequence).join(Genome).join(Taxa) \
                    .order_by(Taxa.genus, Taxa.species, DnaSequence.acc, Bgc.cluster_number) \
                    .all()

    clusters = []
    types = {}

    for entry in ret:
        types[entry.term] = entry.description
        clusters.append({
            "id": "{}_c{}_{}".format(entry.acc, entry.cluster_number, entry.term),
            "parent": entry.term,
            "text": "{} {} Cluster {}".format(entry.species, entry.acc, entry.cluster_number),
            "type": "cluster",
        })

    for name, desc in types.iteritems():
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


@app.route('/api/v1.0/tree/taxa')
def get_taxon_tree():
    '''Get the jsTree structure for all taxa'''
    tree_id = request.args.get('id', '1')

    cur = get_db().cursor()

    HANDLERS = {
        'superkingdom': get_phylum,
        'phylum': get_class,
        'class': get_order,
        'order': get_family,
        'family': get_genus,
        'genus': get_species,
    }

    if tree_id == '1':
        tree = get_superkingdom(cur)
    else:
        params = tree_id.split('_')
        taxlevel = params[0]
        params = params[1:]
        handler = HANDLERS.get(taxlevel, lambda x, y: [])
        tree = handler(cur, params)

    return jsonify(tree)


def get_superkingdom(cur):
    '''Get list of superkingdoms'''
    tree = []
    cur.execute(sql.TAXTREE_SUPERKINGOM)
    kingdoms = cur.fetchall()
    for kingdom in kingdoms:
        tree.append(_create_tree_node('superkingdom_{}'.format(kingdom[0].lower()),
                                      '#', kingdom[0]))

    return tree


def get_phylum(cur, params):
    '''Get list of phyla per kingdom'''
    tree = []
    cur.execute(sql.TAXTREE_PHYLUM, params)
    phyla = cur.fetchall()
    for phylum in phyla:
        id_list = params + [phylum[0].lower()]
        tree.append(_create_tree_node('phylum_{}'.format('_'.join(id_list)),
                                      'superkingdom_{}'.format('_'.join(params)),
                                      phylum[0]))

    return tree


def get_class(cur, params):
    '''Get list of classes per kingdom/phylum'''
    tree = []
    cur.execute(sql.TAXTREE_CLASS, params)
    classes = cur.fetchall()
    for cls in classes:
        id_list = params + [cls[0].lower()]
        tree.append(_create_tree_node('class_{}'.format('_'.join(id_list)),
                                      'phylum_{}'.format('_'.join(params)),
                                      cls[0]))
    return tree


def get_order(cur, params):
    '''Get list of oders per kingdom/phylum/class'''
    tree = []
    cur.execute(sql.TAXTREE_ORDER, params)
    orders = cur.fetchall()
    for order in orders:
        id_list = params + [order[0].lower()]
        tree.append(_create_tree_node('order_{}'.format('_'.join(id_list)),
                                      'class_{}'.format('_'.join(params)),
                                      order[0]))
    return tree


def get_family(cur, params):
    '''Get list of families per kingdom/phylum/class/order'''
    tree = []
    cur.execute(sql.TAXTREE_FAMILY, params)
    families = cur.fetchall()
    for family in families:
        id_list = params + [family[0].lower()]
        tree.append(_create_tree_node('family_{}'.format('_'.join(id_list)),
                                      'order_{}'.format('_'.join(params)),
                                      family[0]))
    return tree


def get_genus(cur, params):
    '''Get list of genera per kingdom/phylum/class/order/family'''
    tree = []
    cur.execute(sql.TAXTREE_GENUS, params)
    genera = cur.fetchall()
    for genus in genera:
        id_list = params + [genus[0].lower()]
        tree.append(_create_tree_node('genus_{}'.format('_'.join(id_list)),
                                      'family_{}'.format('_'.join(params)),
                                      genus[0]))
    return tree


def get_species(cur, params):
    '''Get list of species per kingdom/phylum/class/order/family/genus'''
    tree = []
    cur.execute(sql.TAXTREE_SPECIES, params)
    strains = cur.fetchall()
    for strain in strains:
        tree.append(_create_tree_node('{}'.format(strain.acc.lower()),
                                      'genus_{}'.format('_'.join(params)),
                                      '{} {}.{}'.format(strain.species, strain.acc, strain.version),
                                      disabled=False, leaf=True))

    return tree


@app.route('/api/v1.0/search', methods=['POST'])
def search():
    '''Handle searching the database'''
    search_string = request.json.get('search_string', '')
    try:
        offset = int(request.json.get('offset'))
    except TypeError:
        offset = 0

    try:
        paginate = int(request.json.get('paginate'))
    except TypeError:
        paginate = 50

    total_count, stats, found_bgcs = search_bgcs(search_string, offset=offset, paginate=paginate)

    result = {
        'total': total_count,
        'clusters': found_bgcs,
        'offset': offset,
        'paginate': paginate,
        'stats': stats,
    }

    return jsonify(result)


@app.route('/api/v1.0/export', methods=['POST'])
def export():
    '''Export the search results as CSV file'''
    search_string = request.json.get('search_string', '')
    _, _, found_bgcs = search_bgcs(search_string, mapfunc=create_cluster_csv)

    found_bgcs.insert(0, '#Species\tNCBI accession\tCluster number\tBGC type\tFrom\tTo\tMost similar known cluster\tSimilarity in %\tMIBiG BGC-ID\tResults URL')

    handle = StringIO.StringIO()
    for line in found_bgcs:
        handle.write('{}\n'.format(line))

    handle.seek(0)
    return send_file(handle, attachment_filename='asdb_search_results.csv', as_attachment=True)


@app.route('/api/v1.0/genome/<identifier>')
def show_genome(identifier):
    '''show information for a genome by identifier'''
    search_string = '[acc]{}'.format(identifier)
    _, _, found_bgcs = search_bgcs(search_string)

    return jsonify(found_bgcs)


@app.route('/api/v1.0/available/<category>/<term>')
def list_available(category, term):
    '''list available terms for a given category'''
    return jsonify(available_term_by_category(category, term))


def _create_tree_node(node_id, parent, text, disabled=True, leaf=False):
    '''create a jsTree node structure'''
    ret = {}
    ret['id'] = node_id
    ret['parent'] = parent
    ret['text'] = text
    if disabled:
        ret['state'] = {'disabled': True}
    if leaf:
        ret['type'] = 'strain'
    else:
        ret['children'] = True
    return ret
