import pandas as pd
import requests
from urllib.parse import urlencode
import xml.etree.ElementTree as ET
import time

# ------------------------
# Get result ids from query
# ------------------------

def search_ncbi_proteins(query):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    params = {
        'db': 'protein',
        'term': query,
        'retmode': 'xml',
        'retmax': 100000
    }

    search_url = f"{base_url}?{urlencode(params)}"

    response = requests.get(search_url)
    response.raise_for_status()
    xml_result = response.text

    root = ET.fromstring(xml_result)
    count = root.find('Count').text

    print(f"Total entries found: {count}")
    num_batches = (int(count) + 24) // 25
    print(f"Total available pages: {num_batches}")
    print("(25 entries per page)")

    return xml_result, num_batches

# ------------------------
# Get information of query results (BATCHED version)
# ------------------------

def fetch_protein_info_batch(id_list, batch_number=1):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    BATCH_SIZE = 25

    start_index = (batch_number - 1) * BATCH_SIZE
    end_index = batch_number * BATCH_SIZE
    batch_ids = id_list[start_index:end_index]
    id_string = ",".join(batch_ids)

    params = {
        'db': 'protein',
        'id': id_string,
        'retmode': 'xml'
    }

    fetch_url = f"{base_url}?{urlencode(params)}"

    response = requests.get(fetch_url)
    response.raise_for_status()
    xml_content = response.text

    data = []

    root = ET.fromstring(xml_content)

    def get_text(element, tag_name):
        el = element.find(tag_name)
        return el.text if el is not None else None

    for gbseq in root.findall('.//GBSeq'):
        current_seq = {}

        current_seq['definition'] = get_text(gbseq, 'GBSeq_definition')
        current_seq['primary_accession'] = get_text(gbseq, 'GBSeq_primary-accession')
        current_seq['organism'] = get_text(gbseq, 'GBSeq_organism')
        current_seq['locus'] = get_text(gbseq, 'GBSeq_locus')
        current_seq['length'] = get_text(gbseq, 'GBSeq_length')
        current_seq['moltype'] = get_text(gbseq, 'GBSeq_moltype')
        current_seq['topology'] = get_text(gbseq, 'GBSeq_topology')
        current_seq['division'] = get_text(gbseq, 'GBSeq_division')
        current_seq['update_date'] = get_text(gbseq, 'GBSeq_update-date')
        current_seq['create_date'] = get_text(gbseq, 'GBSeq_create-date')
        current_seq['accession_version'] = get_text(gbseq, 'GBSeq_accession-version')
        current_seq['source'] = get_text(gbseq, 'GBSeq_source')
        current_seq['taxonomy'] = get_text(gbseq, 'GBSeq_taxonomy')

        data.append(current_seq)

    df = pd.DataFrame(data)
    return df

# ------------------------
# Get information from individual entry
# ------------------------

def fetch_single_protein_info(protein_id, export=False, export_folder='files/'):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    params = {
        'db': 'protein',
        'id': protein_id,
        'retmode': 'xml'
    }

    fetch_url = f"{base_url}?{urlencode(params)}"

    response = requests.get(fetch_url)
    response.raise_for_status()

    if export:
        if not export_folder.endswith("/"):
            export_folder = f"{export_folder}/"
        with open(f"{export_folder}{protein_id.replace(".", "_")}_result.xml", "w") as f:
            f.write(response.text)

    return response.text

# ------------------------
# Fetch all data with batching
# ------------------------

def fetch_all_protein_data(id_list, start_batch=1, end_batch=None):
    all_data = []
    
    if end_batch is None:
        end_batch = (len(id_list) + 24) // 25

    for batch_num in range(start_batch, end_batch + 1):
        print(f"Fetching batch {batch_num} of {end_batch}...")
        try:
            batch_df = fetch_protein_info_batch(id_list, batch_num)
            all_data.append(batch_df)
            # Respectful to NCBI servers
            time.sleep(1)
        except Exception as e:
            print(f"Error fetching batch {batch_num}: {e}")
            continue
    
    if all_data:
        return pd.concat(all_data, ignore_index=True)
    else:
        return pd.DataFrame()

# ------------------------
# Test functions
# ------------------------

if __name__ == "__main__":
    query = "Bacillus subtilis[Organism] AND dipicolinate synthase AND subunit A"
    
    # Search NCBI
    protein_results, num_batches = search_ncbi_proteins(query)

    # Get all IDs
    root = ET.fromstring(protein_results)
    id_list = [id_element.text for id_element in root.findall('.//Id')]

    # Get DataFrame for a specific batch range
    start_batch = 1
    end_batch = 3

    # Only fetch a few batches to start with
    df = fetch_all_protein_data(id_list, start_batch, min(end_batch, num_batches))
    
    if not df.empty:
        print(f"Retrieved {len(df)} protein records")
        print(df.head())
        #df.to_csv("files/protein_data.csv", index=False)
    else:
        print("No data retrieved")

    # Get info from individual entry
    protein_id = "SPY12701.1"
    
    # Fetch single protein info
    data = fetch_single_protein_info(protein_id)
    print(f"Retrieved data for {protein_id}")

# END