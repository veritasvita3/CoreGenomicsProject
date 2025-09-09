import requests
import xml.etree.ElementTree as ET
import time
import csv
from collections import defaultdict

def fetch_pubmed_articles(query, retmax=500):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    search_url = f"{base_url}esearch.fcgi"
    params = {
        'db': 'pubmed',
        'term': query,
        'retmax': retmax,
        'retmode': 'json'
    }
    print(f"[DEBUG] Fetching PubMed IDs with query: {query}")
    response = requests.get(search_url, params=params)
    print(f"[DEBUG] PubMed API status code: {response.status_code}")
    data = response.json()
    return data['esearchresult']['idlist']

def fetch_article_details(pubmed_ids):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    fetch_url = f"{base_url}efetch.fcgi"
    
    # Process in batches to avoid URL length issues
    batch_size = 100
    all_articles = []
    
    for i in range(0, len(pubmed_ids), batch_size):
        batch = pubmed_ids[i:i+batch_size]
        params = {
            'db': 'pubmed',
            'id': ','.join(batch),
            'retmode': 'xml'
        }
        print(f"[DEBUG] Fetching details for batch {i//batch_size + 1}")
        response = requests.get(fetch_url, params=params)
        print(f"[DEBUG] efetch status code: {response.status_code}")
        
        # Parse this batch
        articles = parse_articles(response.content)
        all_articles.extend(articles)
        
        # Add delay between batches to respect rate limits
        time.sleep(1)
    
    return all_articles

def parse_articles(xml_data):
    print("[DEBUG] Parsing XML data")
    try:
        root = ET.fromstring(xml_data)
    except ET.ParseError as e:
        print(f"[ERROR] XML parsing failed: {e}")
        return []
    
    articles = []
    for article in root.findall('.//PubmedArticle'):
        # Extract PMID
        pmid_elem = article.find('.//PMID')
        pmid = pmid_elem.text if pmid_elem is not None else "N/A"
        
        # Extract Article Title
        article_title_elem = article.find('.//ArticleTitle')
        article_title = article_title_elem.text if article_title_elem is not None else ""
        
        # Extract DOI
        doi = "N/A"
        article_id_list = article.findall('.//ArticleId')
        for aid in article_id_list:
            if aid.attrib.get('IdType') == 'doi':
                doi = aid.text
                break
        
        # Extract Abstract
        abstract_texts = []
        abstract_elems = article.findall('.//AbstractText')
        for elem in abstract_elems:
            if elem.text is not None:
                abstract_texts.append(elem.text)
        abstract = " ".join(abstract_texts)
        
        # Extract MeSH terms
        mesh_terms = []
        mesh_list = article.findall('.//MeshHeading')
        for mesh in mesh_list:
            descriptor = mesh.find('DescriptorName')
            if descriptor is not None and descriptor.text is not None:
                mesh_terms.append(descriptor.text)
        
        articles.append({
            'pmid': pmid,
            'title': article_title,
            'abstract': abstract,
            'doi': doi,
            'mesh_terms': mesh_terms
        })
    
    print(f"[DEBUG] Parsed {len(articles)} articles")
    return articles

def extract_protein_cancer_associations(articles):
    print(f"[DEBUG] Extracting protein-cancer associations from {len(articles)} articles")
    
    # Dictionary to store proteins and their associated cancers, PMIDs, and DOIs
    protein_data = defaultdict(lambda: {
        'cancers': set(),
        'pmids': set(),
        'dois': set()
    })
    
    # Expanded list of protein-related keywords
    protein_keywords = ['protein', 'antigen', 'receptor', 'enzyme', 'kinase', 'factor', 'peptide', 
                       'ligand', 'antibody', 'cytokine', 'chemokine', 'growth factor', 'hormone',
                       'transcription factor', 'oncogene', 'tumor suppressor']
    
    # Expanded list of cancer types
    cancer_types = [
        'breast cancer', 'lung cancer', 'prostate cancer', 'colorectal cancer', 'leukemia', 
        'melanoma', 'pancreatic cancer', 'ovarian cancer', 'carcinoma', 'sarcoma', 'glioblastoma',
        'bladder cancer', 'liver cancer', 'thyroid cancer', 'cervical cancer', 'endometrial cancer',
        'gastric cancer', 'esophageal cancer', 'head and neck cancer', 'brain cancer', 'lymphoma',
        'myeloma', 'renal cell carcinoma', 'testicular cancer', 'bone cancer', 'skin cancer',
        'oral cancer', 'pharyngeal cancer', 'laryngeal cancer', 'nasopharyngeal cancer'
    ]
    
    for article in articles:
        proteins_found = set()
        cancers_found = set()
        
        # Safely combine title and abstract for text mining
        title = article['title'] or ""
        abstract = article['abstract'] or ""
        text_content = (title + " " + abstract).lower()
        
        # Check MeSH terms for proteins and cancers
        for term in article['mesh_terms']:
            lower_term = term.lower()
            
            # Check for protein mentions in MeSH terms
            if any(keyword in lower_term for keyword in protein_keywords):
                proteins_found.add(term)
            
            # Check for cancer types in MeSH terms
            for cancer in cancer_types:
                if cancer in lower_term:
                    cancers_found.add(cancer)
        
        # Check title and abstract for additional mentions
        for cancer in cancer_types:
            if cancer in text_content:
                cancers_found.add(cancer)
                
        # Look for protein mentions in the text content
        for keyword in protein_keywords:
            if keyword in text_content:
                # Try to extract the protein name
                start_idx = text_content.find(keyword)
                if start_idx != -1:
                    # Extract a window of text around the keyword
                    window_start = max(0, start_idx - 30)
                    window_end = min(len(text_content), start_idx + len(keyword) + 30)
                    text_window = text_content[window_start:window_end]
                    
                    # Add the text window as a potential protein name
                    proteins_found.add(text_window.strip())
        
        # Add to our protein data dictionary
        for protein in proteins_found:
            for cancer in cancers_found:
                protein_data[protein]['cancers'].add(cancer)
                protein_data[protein]['pmids'].add(article['pmid'])
                protein_data[protein]['dois'].add(article['doi'])
    
    # Convert sets to strings for CSV output
    results = []
    for protein, data in protein_data.items():
        if data['cancers']:  # Only include proteins with cancer associations
            results.append({
                'protein': protein,
                'cancer_types': '; '.join(sorted(data['cancers'])),
                'pmids': '; '.join(sorted(data['pmids'])),
                'dois': '; '.join(sorted(data['dois']))
            })
    
    print(f"[DEBUG] Found {len(results)} proteins with cancer associations")
    return results

def main():
    # More specific query to get relevant articles
    query = "(cancer AND (protein OR receptor OR antigen OR enzyme)) AND (breast OR lung OR prostate OR colorectal OR leukemia OR pancreatic OR ovarian OR melanoma)"
    pubmed_ids = fetch_pubmed_articles(query, retmax=500)
    print(f"Found {len(pubmed_ids)} articles")
    
    articles = fetch_article_details(pubmed_ids)
    protein_associations = extract_protein_cancer_associations(articles)
    
    # Sort by number of cancer types (descending)
    protein_associations.sort(key=lambda x: len(x['cancer_types'].split('; ')), reverse=True)
    
    # Limit to 500 unique proteins
    output_associations = protein_associations[:500]
    print(f"[DEBUG] Writing {len(output_associations)} unique proteins to CSV")
    
    # Output to CSV
    with open('cancer_proteins.csv', 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['protein', 'cancer_types', 'pmids', 'dois']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for assoc in output_associations:
            writer.writerow(assoc)
    
    print("Data extraction complete. Output saved to cancer_proteins.csv")

if __name__ == "__main__":
    main()