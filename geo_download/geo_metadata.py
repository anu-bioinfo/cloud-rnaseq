import json
import time
import os
import sys

from Bio import Entrez
import redis
import redis
import hashlib
import re
from ftplib import FTP
import codecs
import argparse

Entrez.email = "yunfang@chanzuckerberg.com"
REDIS_STORE = redis.StrictRedis(host='localhost', port=6379, db=2)


SR_PREFIX = 'sr:pg:'
GDS_PREFIX = 'gds:'
SRA_PREFIX = 'sra:'
PAGE_SIZE = 500
ITEMS_PER_DOWNLOAD = 10



def query_to_sig(query):
    return hashlib.md5(query).hexdigest()

def get_search_results_from_gds(queries):
    '''Getting list of GDS that matches our search query'''
    output_summary = {}
    for query in queries:
        sig = query_to_sig(query)
        page = 0
        page_size = PAGE_SIZE
        while page_size >= PAGE_SIZE:
            try:
                print ("Downloading page %d for query '%s'" % (page, query))
                handle = Entrez.esearch(db="gds", retmax=PAGE_SIZE,
                    retstart=page*PAGE_SIZE, term=query)
                record = Entrez.read(handle)
                handle.close()
                json_str = json.dumps(record['IdList'])
                redis_key = SR_PREFIX + sig + ':' + str(page)
                REDIS_STORE.set(redis_key, json_str)
                page_size = len(record['IdList'])
                page+=1
                time.sleep(1)
            except Exception:
                print("Error getting page %d. Retry in 5 minutes ...." % page)
                time.sleep(5)
        output_summary[sig] = page

def fetch_esummary(id_str, mydb='gds', pfx = GDS_PREFIX):
    while True:
        try:
            handle = Entrez.esummary(db=mydb, id=id_str)
            records = Entrez.read(handle)
            handle.close()

            for r in records:
                key = "%s%s" %(pfx, r['Id'])
                val = json.dumps(r)
                REDIS_STORE.set(key, val)
            return
        except Exception:
            print("ERROR getting %s" % id_str)
            time.sleep(5)

def download_search_results(query, force_download = False):
    query_signature = query_to_sig(query)
    page = 0
    while 1:
        print "Downloading page %d for %s" % (page, query_signature)
        redis_key = SR_PREFIX + query_signature + ':' + str(page)
        res = REDIS_STORE.get(redis_key)
        if res is None:
            return
        id_list = json.loads(res)
        batch_list = []
        for doc_id in id_list:
            if len(batch_list) >= ITEMS_PER_DOWNLOAD:
                #print ("  idx: %d key: %s" % (idx, id_str))
                fetch_esummary(",".join(batch_list), 'gds', GDS_PREFIX)
                time.sleep(1)
                batch_list = []
            if force_download or REDIS_STORE.get(GDS_PREFIX+doc_id) is None:
                batch_list.append(doc_id)
        #print ("  idx: %d key: %s" % (idx, id_str))
        if len(batch_list) > 0:
            fetch_esummary(",".join(batch_list), 'gds', GDS_PREFIX)
        page += 1


sc_re = re.compile("single[\s\-]*cell", re.M | re.U)
def doc_filter(doc_hash):
    if not doc_hash['ExtRelations']:
        return False

    summary = doc_hash['summary'].lower()
    title = doc_hash['title'].lower()
    if not sc_re.search(summary) and not sc_re.search(title):
        return False

    if not doc_hash['PubMedIds']:
        return False

    return True


def get_final_sra_list(queries):
    final_list = []
    for query in queries:
        sig = query_to_sig(query)
        page = 0
        while True:
            redis_key = SR_PREFIX + sig + ':' + str(page)
            res = REDIS_STORE.get(redis_key)
            if res is None:
                break
            id_list = json.loads(res)
            for doc_id in id_list:
                doc_hash = json.loads(REDIS_STORE.get(GDS_PREFIX+doc_id))
                if doc_filter(doc_hash):
                    final_list.append(doc_id)
            page += 1
    return list(set(final_list))

def get_file_list_from_ftp(starting_dir, ftp):
    ftp.cwd(starting_dir)
    flist = ftp.nlst("*/*")
    results = []
    for f in flist:
        results.append(['-', starting_dir + '/' + f, 0])
    return results

def get_sra_files(doc_id):
    redis_key = GDS_PREFIX + doc_id
    rec = json.loads(REDIS_STORE.get(redis_key))
    sra_list = rec['ExtRelations']
    result_list = []
    idx = 0
    tries = 0
    while idx < len(sra_list):
        s = sra_list[idx]
        if s["RelationType"] != 'SRA':
            continue
        sra_url = s["TargetFTPLink"]
        m = re.match("ftp\:\/\/([^\/]+)(\/.*)", sra_url)
        try:
            ftp_host = m.group(1)
            ftp_dir  = m.group(2)
            ftph = FTP(ftp_host)
            ftph.login()
            result_list = result_list + get_file_list_from_ftp(ftp_dir, ftph)
            ftph.close()
            tries = 0
            idx += 1
        except Exception:
            tries += 1
            if tries < 3:
                print("ERROR getting %s. Retry in 5 seconds..." % sra_url)
                time.sleep(5)
            else:
                idx += 1 # failed too many times. move on to the next one.
                tries = 0

    return result_list

def download_sra_urls(doc_list, force_download = False):
    for doc_id in doc_list:
        print("=======Downloading SRA for doc %s ====" % doc_id)
        redis_key = SRA_PREFIX + doc_id
        if force_download or REDIS_STORE.get(redis_key) is None:
            res = get_sra_files(doc_id)
            redis_val = json.dumps(res)
            REDIS_STORE.set(redis_key, redis_val)

def doc_list_to_file(doc_list, filename):
    with open(filename, 'w') as fd:
        for doc_id in doc_list:
            fd.write("%s\n" % doc_id)
def doc_summary_to_file(doc_list, filename):
    with codecs.open(filename,'w', encoding='utf-8') as fd:
        idx = 1
        for doc_id in doc_list:
            redis_key = 'gds:' + doc_id
            rec = json.loads(REDIS_STORE.get(redis_key))
            pubmed_str = ",".join([str(x) for x in rec['PubMedIds']])
            output = "%d. Doc ID: %s, GSE: %s: PubMed Id: %s, Taxon: %s, Date: %s, Samples: %s\n" % \
                (idx, doc_id, str(rec['GSE']), pubmed_str,
                 rec['taxon'], rec['PDAT'], str(rec['n_samples']))

            output += "Title: %s \n" % rec['title']
            output += "Summary: %s \n" % rec['summary']
            output += "# SRA URLS: %d \n" % len(rec['ExtRelations'])
            output += "=============================================================\n"
            fd.write("%s\n" % output)
            idx += 1
def main():
    parser = argparse.ArgumentParser(
        description='Download metadata from GEO datasets')
    parser.add_argument('-f', action="store", dest='query_list_file', default=False)
    results = parser.parse_args()
    if results.query_list_file:
        queries = []
        with open(results.query_list_file, 'r') as query_list_f:
            for line in query_list_f:
                queries.append(line.rstrip())

        get_search_results_from_gds(queries)

        for query in queries:
            download_search_results(query)

        final_download_list = get_final_sra_list(queries)
        download_sra_urls(final_download_list)
        doc_list_to_file(final_download_list, 'sra_doc_ids.txt')
        doc_summary_to_file(final_download_list, 'sra_doc_summaries.txt')
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
