"""Minimal stdlib RCSB / UniProt client with on-disk JSON cache (single-threaded, polite)."""
import hashlib
import json
import os
import time
import urllib.error
import urllib.request

UA = {'User-Agent': 'laca-imprint-audit/1.0 (academic)', 'Content-Type': 'application/json'}
CACHE = None


def set_cache(d):
    global CACHE
    CACHE = d
    os.makedirs(d, exist_ok=True)


def _cached(key, fn):
    h = hashlib.sha1(key.encode()).hexdigest()
    p = os.path.join(CACHE, h + '.json')
    if os.path.exists(p):
        return json.load(open(p))
    out = fn()
    tmp = p + '.tmp'
    json.dump(out, open(tmp, 'w'))
    os.replace(tmp, p)
    return out


def _req(url, data=None, tries=5):
    for k in range(tries):
        try:
            r = urllib.request.Request(url, data=data, headers=UA)
            with urllib.request.urlopen(r, timeout=120) as f:
                raw = f.read()
            time.sleep(0.15)
            return raw
        except urllib.error.HTTPError as e:
            if e.code == 204:
                return b''
            if e.code in (429, 500, 502, 503, 504) and k < tries - 1:
                time.sleep(2 * (k + 1)); continue
            raise
        except Exception:
            if k < tries - 1:
                time.sleep(2 * (k + 1)); continue
            raise


def search(query, return_type='entry'):
    body = {'query': query, 'return_type': return_type, 'request_options': {'return_all_hits': True}}
    key = 'search:' + json.dumps(body, sort_keys=True)

    def fn():
        raw = _req('https://search.rcsb.org/rcsbsearch/v2/query', json.dumps(body).encode())
        if not raw:
            return []
        return [x['identifier'] for x in json.loads(raw).get('result_set', [])]
    return _cached(key, fn)


def graphql(q):
    key = 'gql:' + q

    def fn():
        raw = _req('https://data.rcsb.org/graphql', json.dumps({'query': q}).encode())
        return json.loads(raw)
    return _cached(key, fn)


def uniprot(acc):
    key = 'uniprot:' + acc

    def fn():
        try:
            raw = _req(f'https://rest.uniprot.org/uniprotkb/{acc}.json')
            return json.loads(raw)
        except urllib.error.HTTPError as e:
            return {'_error': e.code}
    return _cached(key, fn)


def download_cif(pdb, outdir):
    p = os.path.join(outdir, pdb.upper() + '.cif')
    if os.path.exists(p) and os.path.getsize(p) > 0:
        return p
    raw = _req(f'https://files.rcsb.org/download/{pdb.upper()}.cif')
    tmp = p + '.tmp'
    open(tmp, 'wb').write(raw)
    os.replace(tmp, p)
    return p


def terminal(attr, op, value, service='text'):
    return {'type': 'terminal', 'service': service, 'parameters': {'attribute': attr, 'operator': op, 'value': value}}


def AND(*nodes):
    return {'type': 'group', 'logical_operator': 'and', 'nodes': list(nodes)}
