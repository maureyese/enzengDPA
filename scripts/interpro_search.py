import requests
import json
import time

def run_interpro(email: str, fasta: str) -> str:
    """
    Send input to InterProScan and get job ID.
    
    Args:
        email: Email for notification
        fasta: FASTA sequence
    
    Returns:
        Job ID as string
    """
    search_url = 'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run'
    
    # Prepare parameters
    search_params = {
        'email': email,
        'sequence': fasta
    }
    
    try:
        # Perform POST request
        response = requests.post(
            search_url,
            headers={
                "Content-Type": "application/x-www-form-urlencoded",
                "Accept": "text/plain"
            },
            data=search_params
        )
        
        # Validate response
        response.raise_for_status()
        
        # Get run ID
        job_id = response.text.strip()
        return job_id
        
    except requests.exceptions.RequestException as error:
        print(f"Error running InterProScan: {error}")
        raise error

def check_interpro_status(job_id: str) -> str:
    """
    Check status of InterProScan job.
    
    Args:
        job_id: Job ID to check
    
    Returns:
        Job status as string
    """
    status_url = f'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/{job_id}'
    
    try:
        response = requests.get(
            status_url,
            headers={"Accept": "text/plain"}
        )
        
        # Validate response
        response.raise_for_status()
        
        # Get job status
        job_status = response.text.strip()
        return job_status
        
    except requests.exceptions.RequestException as error:
        print(f"Error checking job status: {error}")
        raise error

def get_interpro_results(job_id: str, output_file: str = None) -> dict:
    """
    Retrieve results from InterProScan job.
    
    Args:
        job_id: Job ID to retrieve results for
        output_file: Optional filename to save JSON results
    
    Returns:
        Results as dictionary
    """
    results_url = f'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{job_id}/json'
    
    try:
        # Request JSON results
        response = requests.get(
            results_url,
            headers={"Accept": "application/json"}
        )
        
        # Validate response
        response.raise_for_status()
        
        # Get JSON results
        results = response.json()
        
        # Export file if output_file is specified
        if output_file:
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(results, f, indent=2)
            print(f"Results saved to {output_file}")
        
        return results
        
    except requests.exceptions.RequestException as error:
        print(f"Error fetching InterProScan results: {error}")
        raise error
    except json.JSONDecodeError as error:
        print(f"Error parsing JSON response: {error}")
        raise error

# ------------------------
# Testing functions
# ------------------------

if __name__ == "__main__":
    fasta = """>AAA41145.1 fatty acid synthase [Rattus norvegicus]
MEEVVIAGMSGKLPESENLQEFWANLIGGVDMVTDDDRRWKAGLYGLPKRSGKLKDLSKFDASFFGVHPK
QAHTMDPQLRLLLEVSYEAIVDGGINPASLRGTNTGVWVGVSGSEASEALSRDPETLLGYSMVGCQRAMM
ANRLSFFFDFKGPSIALDTACSSSLLALQNAYQAIRSGECPAATVGGINLLLKPNTSVQFMKLGMLSPDG
TCRSFDDSGNGYCRAEAVVAVLLTKKSLARRVYATILNAGTNTDGCKEQGVTFPSGEAQEQLIRSLYQPG
GVAPESLEYIEAHGTGTKVGDPQELNGITRSLCAFRQSPLLIGSTKSNMGHPEPASGLAALTKVLLSLEN
GVWAPNLHFHNPNPEIPALLDGRLQVVDRPLPVRGGIVGINSFGFGGANVHVILQPNTQQAPAPAPHAAL
PHLLHASGRTMEAVQGLLEQGRQHSQDLAFVSMLNDIAATPTAAMPFRGYTVLGVEGHVQEVQQVPASQR
PLWFICSGMGTQWRGMGLSLMRLDSFRESILRSDEALKPLGVKVSDLLLSTDEHTFDDIVHSFVSLTAIQ
IALIDLLTSMGLKPDGIIGHSLGEVACGYADGCLSQREAVLAAYWRGQCIKDANLPAGSMAAVGLSWEEC
KQRCPPGVVPACHNSEDTVTISGPQAAVNEFVEQLKQEGVFAKEVRTGGLAFHSYFMEGIAPTLLQALKK
VIREPRPRSARWLSTSIPEAQWQSSLARTSSAEYNVNNLVSPVLFQEALWHVPEHAVVLEIAPHALLQAV
LKRGVKPSCTIIPLMKRDHKDNLEFFLTNLGKVHLTGIDINPNALFPPVEFPVPRGTPLISPHIKWDHSQ
TWDIPVAEDFPNGSSSSSATVYNIDASSESSDHYLVDHCIDGRVLFPGTGYLYLVWKTLARSLSLSLEET
PVVFENVTFHQATILPRTGTVPLEVRLLEASHAFEVSDSGNLIVSGKVYQWEDPDSKLFDHPEVPIPAES
ESVSRLTQGEVYKELRLRGYDYGPHFQGVYEATLEGEQGKLLWKDNWVTFMDTMLQISILGFSKQSLQLP
TRVTAIYIDPATHLQKVYMLEGDTQVADVTTSRCLGVTVSGGVYISRLQTTATSRRQQEQLVPTLEKFVF
TPHVEPECLSESAILQKELQLCKGLAKALQTKATQQGLKMTVPGLEDLPQHGLPRLLAAACQLQLNGNLQ
LELGEVLARERLLLPEDPLISGLLNSQALKACIDTALENLSTLKMKVVEVLAGEGHLYSHISALLNTQPM
LQLEYTATDRHPQALKDVQTKLQQHDVAQGQWDPSGPAPTNLGALDLVVCNCALATLGDPALALDNMVAA
LKDGGFLLMHTVLKGHALGETLACLPSEVQPGPSFLSQEEWESLFSRKALHLVGLKKSFYGTALFLCRRL
SPQDKPIFLPVEDTSFQWVDSLKSILATSSSQPVWLTAMNCPTSGVVGLVNCLRKEPGGHRIRCILLSNL
SSTSHVPKLDPGSSELQKVLESDLVMNVYRDGAWGAFRHFQLEQDKPEEQTAHAFVNVLTRGDLASIRWV
SSPLKHMQPPSSSGAQLCTVYYASLNFRDIMLATGKLSPDAIPGKWASRDCMLGMEFSGRDKCGRRVMGL
VPAEGLATSVLLSPDFLWDVPSSWTLEEAASVPVVYTTAYYSLVVRGRIQHGETVLIHSGSGGVGQAAIS
IALSLGCRVFTTVGSAEKRAYLQARFPQLDDTSFANSRDTSFEQHVLLHTGGKGVDLVLNSLAEEKLQAS
VRCLAQHGRFLEIGKFDLSNNHPLGMAIFLKNVTFHGILLDALFEGANDSWREVAELLKAGIRDGVVKPL
KCTVFPKAQVEDAFRYMAQGKHIGKVLVQVREEEPEAMLPGAQPTLISAISKTFCPEHKSYIITGGLGGF
GLELARWLVLRGAQRLVLTSRSGIRTGYQAKHVREWRRQGIHVLVSTSNVSSLEGARALIAEATKLGPVG
GVFNLAMVLRDAMLENQTPELFQDVNKPKYNGTLNLDRATREACPELDYFVAFSSVSCGRGNAGQSNYGF
ANSTMERICEQRRHDGLPGLAVQWGAIGDVGIILEAMGTNDTVVGGTLPQRISSCMEVLDLFLNQPHAVL
SSFVLVEKKAVAHGDGEAQRDLVKAVAHILGIRDLAGINLDSSLADLGLDSLMGVEVRQILEREHDLVLP
IREVRQLTLRKLQEMSSKAGSDTELAAPKSKNDTSLKQAQLNLSILLVNPEGPTLTRLNSVQSSERPLFL
VHPIEGSITVFHSLAAKLSVPTYGLQCTQAAPLDSIPNLAAYYIDCIKQVQPEGPYRVAGYSFGACVAFE
MCSQLQAQQGPAPAHNNLFLFDGSHTYVLAYTQSYRAKLTPGCEAEAEAEAICFFIKQFVDAEHSKVLEA
LLPLKSLEDRVAAAVDLITRSHQSLDRRDLSFAAVSFYYKLRAADQYKPKAKYHGNVILLRAKTGGTYGE
DLGADYNLSQVCDGKVSVHIIEGDHRTLLEGRGLESIINIIHSSLAEPRVSVREG"""
    
    email = 'example@example.com'  # E-Mail for any issue
    
    try:
        # Submit job
        job_id = run_interpro(email, fasta)
        print(f"Job ID: {job_id}")
        
        # Poll for status
        while True:
            print(f"Checking job status from {job_id}...")
            job_status = check_interpro_status(job_id)
            print(f"Current Status: {job_status}")
            
            if job_status in ["QUEUED", "RUNNING"]:
                time.sleep(5)  # Wait 5 seconds before checking again
                continue
            elif job_status == "FINISHED":
                print("Job completed successfully!")
                break
            else:
                print(f"Job failed with status: {job_status}")
        
        # Get results
        output_file = "interPro_results.json"
        results = get_interpro_results(job_id, output_file)
        
        print("Results retrieved successfully!")
        # Optional: print a summary of results
        if results and 'results' in results:
            print(f"Found {len(results['results'])} domain matches")
        
    except Exception as error:
        print(f"Error running and checking InterProScan job: {error}")
