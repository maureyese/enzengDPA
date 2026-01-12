import requests
import time
from urllib.parse import urlencode
import pandas as pd
import xml.etree.ElementTree as ET
from io import StringIO

# ------------------------
# Run BLAST
# ------------------------

def run_blast(
    sequence:str,
    organism:str,
    database:str,
    program:str,
    word_size: int = 6,          # Default for blastp is 6
    expect_value: float = 10.0,  # Default for blastp is 10
    matrix: str = 'BLOSUM62',    # Default matrix
    gapcosts: str = '11 1',      # Default gap costs (e.g., 11 for gap open, 1 for gap extend)
    filter_string: str = 'L',    # 'L' for Low-compositional complexity filter (like F in the old format)
    composition_stats: int = 2,  # Default for composition-based statistics (2 or 1)
    export:bool = False,
    export_folder = "files/"
) -> tuple[str | None, str | None]:
    '''Run BLAST through API guidelines'''
    # Set URL for BLAST
    BASE_URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'

    # 1) Send the 'Put' request to start the BLAST job
    params_put = {
        'CMD': 'Put',
        'QUERY': sequence,
        'DATABASE': 'nr',  # 'nr' for blastx, or e.g. 'nt' for blastn
        'PROGRAM': 'blastp', # 'blastp', 'blastn', 'blastx'
        'ENTREZ_QUERY': organism,
        'OUTPUT_TYPE': 'XML',
        # ---- Additional parameters ----
        'WORD_SIZE': word_size,
        'EXPECT': expect_value,
        'MATRIX': matrix,
        'GAPCOSTS': gapcosts,
        'FILTER': filter_string,
        'COMPOSITION_BASED_STATISTICS': composition_stats,
    }

    put_url = f"{BASE_URL}?{urlencode(params_put)}"

    try:
        # Perform the PUT request
        response_put = requests.get(put_url)
        response_put.raise_for_status()
        text_put = response_put.text

        # Extract the RID (Request ID) from the response text
        if 'RID = ' in text_put:
            rid = text_put.split('RID = ')[1].split('\n')[0].strip()
            print('RID obtained:', rid)

            # 2) Poll for the job status in a loop
            is_ready = False
            poll_count = 0
            last_status = ''

            while not is_ready:
                poll_count += 1

                # Build the URL for status check
                params_get = {
                    'CMD': 'Get',
                    'RID': rid,
                    'FORMAT_OBJECT': 'SearchInfo',
                }
                get_url = f"{BASE_URL}?{urlencode(params_get)}"

                response_get = requests.get(get_url)
                response_get.raise_for_status()
                text_get = response_get.text

                if 'Status=WAITING' in text_get:
                    last_status = 'WAITING'
                    print(f"({poll_count}) BLAST is still running. Waiting 15 seconds...")
                    time.sleep(15)

                elif 'Status=FAILED' in text_get:
                    last_status = 'FAILED'
                    print(f"({poll_count}) BLAST job failed.")
                    break

                elif 'Status=UNKNOWN' in text_get:
                    last_status = 'UNKNOWN'
                    print(f"({poll_count}) BLAST job unknown (possibly expired or invalid RID).")
                    break

                elif 'Status=READY' in text_get:
                    last_status = 'READY'

                    if 'ThereAreHits=yes' in text_get:
                        print(f"({poll_count}) BLAST job is complete, and hits are found!")
                    else:
                        print(f"({poll_count}) BLAST job is complete, but NO hits found.")

                    is_ready = True
                    break

                else:
                    last_status = 'OTHER'
                    print(f"({poll_count}) Unexpected status. Stopping...")
                    break

            # 3) If the job completed successfully, retrieve the final XML
            if last_status == 'READY':
                params_result = {
                    'CMD': 'Get',
                    'RID': rid,
                    'FORMAT_TYPE': 'XML',
                }
                result_url = f"{BASE_URL}?{urlencode(params_result)}"

                response_result = requests.get(result_url)
                response_result.raise_for_status()
                xml_result = response_result.text

                print('=== BLAST XML RESULTS ===\n')
                #print(xml_result)

                # Save to file
                if export:
                    if not export_folder.endswith("/"):
                        export_folder = f"{export_folder}/"
                    with open(f"{export_folder}{rid}_results.xml", "w") as f:
                        f.write(xml_result)
                    print(f"Results written to {rid}_results.xml")
                else:
                    print(f"BLAST Results not exported (Argument 'export' = False). RID: {rid}\n")

                return rid, xml_result
            else:
                return None, None
        else:
            raise Exception('No RID found in the PUT response. Aborting.')

    except Exception as error:
        print('Error:', str(error))
        return None, None

# ------------------------
# Get BLAST results using RID
# ------------------------

def get_blast_results(
        rid: str,
        export:bool = False,
        export_folder:str = "files/"
) -> tuple[str | None, str | None]:
    '''Get BLAST results using RID'''
    BASE_URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'

    is_ready = False
    poll_count = 0
    last_status = ''

    while not is_ready:
        poll_count += 1

        params_get = {
            'CMD': 'Get',
            'RID': rid,
            'FORMAT_OBJECT': 'SearchInfo',
        }
        get_url = f"{BASE_URL}?{urlencode(params_get)}"

        response_get = requests.get(get_url)
        response_get.raise_for_status()
        text_get = response_get.text

        if 'Status=WAITING' in text_get:
            last_status = 'WAITING'
            print(f"({poll_count}) BLAST is still running. Waiting 15 seconds...")
            time.sleep(15)

        elif 'Status=FAILED' in text_get:
            last_status = 'FAILED'
            print(f"({poll_count}) BLAST job failed.")
            break

        elif 'Status=UNKNOWN' in text_get:
            last_status = 'UNKNOWN'
            print(f"({poll_count}) BLAST job unknown (possibly expired or invalid RID).")
            break

        elif 'Status=READY' in text_get:
            last_status = 'READY'

            if 'ThereAreHits=yes' in text_get:
                print(f"({poll_count}) BLAST job is complete, and hits are found!")
            else:
                print(f"({poll_count}) BLAST job is complete, but NO hits found.")

            is_ready = True
            break

        else:
            last_status = 'OTHER'
            print(f"({poll_count}) Unexpected status. Stopping...")
            break

    # 3) If the job completed successfully, retrieve the final XML
    if last_status == 'READY':
        params_result = {
            'CMD': 'Get',
            'RID': rid,
            'FORMAT_TYPE': 'XML',
        }
        result_url = f"{BASE_URL}?{urlencode(params_result)}"

        response_result = requests.get(result_url)
        response_result.raise_for_status()
        xml_result = response_result.text

        print('=== BLAST XML RESULTS ===\n')
        #print(xml_result)

        # Save to file
        if export:
            if not export_folder.endswith("/"):
                export_folder = f"{export_folder}/"
            with open(f"{export_folder}{rid}_results.xml", "w") as f:
                f.write(xml_result)
            print(f"Results written to {rid}_results.xml")
        else:
            print(f"BLAST Results not exported (Argument 'export' = False). RID: {rid}")

        return rid, xml_result
    else:
        return None, None

# ------------------------
# Convert BLAST XML fie to DataFrame
# ------------------------

def convert_blast_xml_to_pd(
        xml_results: str, 
        rid: str,
        query_data: str = None,
        export: bool = False,
        export_folder: str = "files/"
) -> pd.DataFrame:
    '''Convert BLAST XML results to Pandas DataFrame and prepend query info'''
    
    # 1. Parse BLAST Results
    root = ET.fromstring(xml_results)
    data = []
    namespace = ''

    for hit in root.findall(f'.//{namespace}Hit'):
        hit_num = hit.find(f'{namespace}Hit_num').text
        hit_id = hit.find(f'{namespace}Hit_id').text
        hit_def = hit.find(f'{namespace}Hit_def').text
        hit_accession = hit.find(f'{namespace}Hit_accession').text
        hit_len = hit.find(f'{namespace}Hit_len').text

        hsp = hit.find(f'.//{namespace}Hsp')
        if hsp is not None:
            identity = int(hsp.find(f'{namespace}Hsp_identity').text)
            align_len = int(hsp.find(f'{namespace}Hsp_align-len').text)
            
            data.append({
                'Hit_num': hit_num,
                'Hit_id': hit_id,
                'Hit_def': hit_def,
                'Hit_accession': hit_accession,
                'Hit_len': hit_len,
                'Hsp_num': hsp.find(f'{namespace}Hsp_num').text,
                'Bit_score': float(hsp.find(f'{namespace}Hsp_bit-score').text),
                'Score': int(hsp.find(f'{namespace}Hsp_score').text),
                'E_value': hsp.find(f'{namespace}Hsp_evalue').text,
                'Query_from': int(hsp.find(f'{namespace}Hsp_query-from').text),
                'Query_to': int(hsp.find(f'{namespace}Hsp_query-to').text),
                'Hit_from': int(hsp.find(f'{namespace}Hsp_hit-from').text),
                'Hit_to': int(hsp.find(f'{namespace}Hsp_hit-to').text),
                'Identity': identity,
                'Positive': int(hsp.find(f'{namespace}Hsp_positive').text),
                'Gaps': int(hsp.find(f'{namespace}Hsp_gaps').text),
                'Align_len': align_len,
                'Percent_identity': round(identity / align_len * 100, 2),
                'Query_seq': hsp.find(f'{namespace}Hsp_qseq').text,
                'Hit_seq': hsp.find(f'{namespace}Hsp_hseq').text,
                'Midline': hsp.find(f'{namespace}Hsp_midline').text,
            })

    df = pd.DataFrame(data)

    # 2. Parse Query Data (GenBank XML) and Prepend
    if query_data:
        q_root = ET.fromstring(query_data)
        # Handle finding the sequence and definition in the GBSet format
        q_seq_elem = q_root.find('.//GBSeq_sequence')
        q_def_elem = q_root.find('.//GBSeq_definition')
        q_acc_elem = q_root.find('.//GBSeq_accession-version')
        
        if q_seq_elem is not None:
            # Create a row representing the query itself
            query_row = pd.DataFrame([{
                'Hit_num': '0',
                'Hit_id': 'Query',
                'Hit_def': q_def_elem.text if q_def_elem is not None else 'Query Sequence',
                'Hit_accession': q_acc_elem.text if q_acc_elem is not None else rid,
                'Hit_len': len(q_seq_elem.text),
                'Bit_score': 0.0,
                'E_value': '0',
                'Percent_identity': 100.0,
                'Query_seq': q_seq_elem.text.upper(),
                'Hit_seq': q_seq_elem.text.upper(),
            }])
            
            # Prepend to the dataframe
            df = pd.concat([query_row, df], ignore_index=True).fillna('')
    if export:
        if not export_folder.endswith("/"):
            export_folder = f"{export_folder}/"
        df.to_csv(f"{export_folder}{rid}_results.csv", index=False)
        print(f"DataFrame saved to {export_folder}{rid}_results.csv")

    return df

# ------------------------
# Export BLAST results to FASTA file
# ------------------------

def export_blast_results_to_fasta(
        rid: str,
        df: pd.DataFrame,
        export_folder: str = "files/"
) -> None:
    '''Create FASTA file from BLASTA DataFrame'''
    try:
        if df.empty:
            print("DataFrame is empty. Nothing to export.")
            return

        if not export_folder.endswith("/"):
            export_folder = f"{export_folder}/"
        
        output_filename = f"{export_folder}{rid}_blast.fasta"

        with open(output_filename, 'w') as fasta_file:
            for index, row in df.iterrows():
                accession = str(row.get('Hit_accession', '')).strip()
                if not accession or accession == 'nan':
                    accession = str(row.get('Hit_id', 'Unknown'))
                
                header = row.get('Hit_def', 'No definition')
                # --- FIX END ---

                if index == 0:
                    sequence = row.get('Query_seq', '')
                else:
                    sequence = row.get('Hit_seq', '')

                if sequence:
                    fasta_file.write(f">{accession} {header}\n")
                    fasta_file.write(f"{sequence}\n")

        print(f"BLAST results successfully exported to {output_filename}")
        return output_filename
        
    except Exception as e:
        print(f"Error exporting BLAST results: {e}")

# ------------------------
# Testing functions
# ------------------------

if __name__ == "__main__":
    sequence_id = "msslkgkrigfgltgshctyeavfpqievlvnegaevrpvvtfnvkstntrfgegaewvkkieeltgyeaidsivkaeplgpklpldcmviapltgnsmsklanamtdspvlmaakatirnnrpvvlgistndalglngtnlmrlmstkniffipfgqddpfkkpnsmvakmdllpqtiekalmhqqlqpilvenyqgnd"
    organism = "Bacillus[Organism]"
    database = "nr_cluster_seq"
    program = "blastp"
    # Perform BLAST
    if sequence_id:
        rid, xml_results = run_blast(
            sequence=sequence_id,
            organism=organism,
            database=database,
            program=program
        )
    else:
        print("Please provide a sequence")
    
    # Try getting XML results from RID
    rid, xml_results = get_blast_results(rid=rid)
    # Convert results to DataFrame
    blast_df = convert_blast_xml_to_pd(xml_results=xml_results, rid=rid, export=True)
    # Export BLAST results as FASTA file
    #export_blast_results_to_fasta(rid=rid, query_id="SPY12701.1", df=blast_df)

# END