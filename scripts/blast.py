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
        xml_results:str, 
        rid:str,
        export:bool = False,
        export_folder:str = "files/"
) -> pd.DataFrame:
    '''Convert BLAST XML results to Pandas DataFrame'''
    # Parse XML
    root = ET.fromstring(xml_results)

    # Extract data
    data = []

    # Namespace for BLAST XML
    namespace = ''

    # Extract data from each hit
    for hit in root.findall(f'.//{namespace}Hit'):
        hit_num = hit.find(f'{namespace}Hit_num').text
        hit_id = hit.find(f'{namespace}Hit_id').text
        hit_def = hit.find(f'{namespace}Hit_def').text
        hit_accession = hit.find(f'{namespace}Hit_accession').text
        hit_len = hit.find(f'{namespace}Hit_len').text

        # Get HSP data
        hsp = hit.find(f'.//{namespace}Hsp')
        if hsp is not None:
            hsp_num = hsp.find(f'{namespace}Hsp_num').text
            bit_score = hsp.find(f'{namespace}Hsp_bit-score').text
            score = hsp.find(f'{namespace}Hsp_score').text
            evalue = hsp.find(f'{namespace}Hsp_evalue').text
            query_from = hsp.find(f'{namespace}Hsp_query-from').text
            query_to = hsp.find(f'{namespace}Hsp_query-to').text
            hit_from = hsp.find(f'{namespace}Hsp_hit-from').text
            hit_to = hsp.find(f'{namespace}Hsp_hit-to').text
            identity = hsp.find(f'{namespace}Hsp_identity').text
            positive = hsp.find(f'{namespace}Hsp_positive').text
            gaps = hsp.find(f'{namespace}Hsp_gaps').text
            align_len = hsp.find(f'{namespace}Hsp_align-len').text
            qseq = hsp.find(f'{namespace}Hsp_qseq').text
            hseq = hsp.find(f'{namespace}Hsp_hseq').text
            midline = hsp.find(f'{namespace}Hsp_midline').text

            # Append to data list
            data.append({
                'Hit_num': hit_num,
                'Hit_id': hit_id,
                'Hit_def': hit_def,
                'Hit_accession': hit_accession,
                'Hit_len': hit_len,
                'Hsp_num': hsp_num,
                'Bit_score': float(bit_score),
                'Score': int(score),
                'E_value': evalue,
                'Query_from': int(query_from),
                'Query_to': int(query_to),
                'Hit_from': int(hit_from),
                'Hit_to': int(hit_to),
                'Identity': int(identity),
                'Positive': int(positive),
                'Gaps': int(gaps),
                'Align_len': int(align_len),
                'Percent_identity': round(int(identity) / int(align_len) * 100, 2),
                'Query_seq': qseq,
                'Hit_seq': hseq,
                'Midline': midline,
            })

    # Create DataFrame
    df = pd.DataFrame(data)

    # Display the DataFrame info and first few rows
    print(f"DataFrame shape: {df.shape}")
    #print("\nFirst 5 rows:")
    #print(df.head())

    # Save dataframe
                    # Save to file
    if export:
        if not export_folder.endswith("/"):
            export_folder = f"{export_folder}/"
        df.to_csv(f"{export_folder}{rid}_results.csv", index=False)
        print(f"DataFrame saved to {rid}_results.csv")

    return df

# ------------------------
# Export BLAST results to FASTA file
# ------------------------

def export_blast_results_to_fasta(
        rid:str,
        query_id:str,
        df: pd.DataFrame,
        export_folder:str = "files/"
) -> None:
    try:
        if not export_folder.endswith("/"):
            export_folder = f"{export_folder}/"
        # Create output filename
        output_filename = f"{export_folder}{rid}_blast.fasta"

        # Open file for writing
        with open(output_filename, 'w') as fasta_file:
            # Get the template sequence
            template = df['Query_seq'].iloc(0)
            fasta_file.write(f">{query_id}\n")
            fasta_file.write(f"{template}\n")

            for index, row in df.iterrows():
                # Get title from Hit_def column
                title = row['Hit_def']
                # Get sequence from Hit_seq
                sequence = row['Hit_seq']

                # Write FASTA entry directly to file
                fasta_file.write(f">{title}\n")
                fasta_file.write(f"{sequence}\n")

        print(f"BLAST results successfully exported to {output_filename}")
        print(f"Total sequences exported: {len(df)}")
    except KeyError as e:
        print(f"Error: Missing required column in DataFrame: {e}")
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