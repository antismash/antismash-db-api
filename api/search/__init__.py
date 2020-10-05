'''Search-related functions'''
from sqlalchemy import (
    func,
)
from api.models import (
    db,
    AsDomain,
    BgcType,
    Cds,
    DnaSequence,
    Genome,
    Region,
    Taxa,
    t_rel_regions_types,
)
from .clusters import (
    cluster_query_from_term,
    CLUSTER_FORMATTERS,
)
from .genes import (
    gene_query_from_term,
    GENE_FORMATTERS,
)
from .domains import (
    domain_query_from_term,
    DOMAIN_FORMATTERS,
)

from .helpers import UnknownQueryError

#######
# The following imports are just so the code depending on search doesn't need changes
from .available import available_term_by_category  # noqa: F401
#######

FORMATTERS = {
    'cluster': CLUSTER_FORMATTERS,
    'gene': GENE_FORMATTERS,
    'domain': DOMAIN_FORMATTERS,
}


class NoneQuery(object):
    '''A 'no result' return object'''
    def all(self):
        '''Just return an empty list'''
        return []


def core_search(query):
    '''Actually run the search logic'''
    sql_query = NoneQuery()

    if query.search_type == 'cluster':
        sql_query = cluster_query_from_term(query.terms).order_by(Region.region_id)
    elif query.search_type == 'gene':
        sql_query = gene_query_from_term(query.terms).order_by(Cds.cds_id)
    elif query.search_type == 'domain':
        sql_query = domain_query_from_term(query.terms).order_by(AsDomain.as_domain_id)
    else:
        raise UnknownQueryError()

    results = sql_query.all()

    return results


def format_results(query, results):
    '''Get the appropriate formatter for the query'''
    try:
        fmt_func = FORMATTERS[query.search_type][query.return_type]
        return fmt_func(results)
    except KeyError:
        return []


def region_stats(regions):
    '''Calculate stats on the search results'''
    stats = {}

    if len(regions) < 1:
        return stats

    bgc_ids = set(map(lambda x: x.region_id, regions))

    clusters_by_type_list = db.session.query(BgcType.term, func.count(BgcType.term)) \
                                      .join(t_rel_regions_types).join(Region) \
                                      .filter(Region.region_id.in_(bgc_ids)).group_by(BgcType.term).order_by(BgcType.term).all()
    clusters_by_type = {}
    if clusters_by_type_list is not None:
        clusters_by_type['labels'], clusters_by_type['data'] = zip(*clusters_by_type_list)
    stats['clusters_by_type'] = clusters_by_type

    clusters_by_phylum_list = db.session.query(Taxa.phylum, func.count(Taxa.phylum)) \
                                        .join(Genome).join(DnaSequence).join(Region) \
                                        .filter(Region.region_id.in_(bgc_ids)).group_by(Taxa.phylum).all()
    clusters_by_phylum = {}
    if clusters_by_phylum_list is not None:
        clusters_by_phylum['labels'], clusters_by_phylum['data'] = zip(*clusters_by_phylum_list)
    stats['clusters_by_phylum'] = clusters_by_phylum

    return stats
