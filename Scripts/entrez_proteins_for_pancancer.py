import requests
import xml.etree.ElementTree as ET
import time
import csv

def fetch_pubmed_articles(query, retmax=300):
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
    params = {
        'db': 'pubmed',
        'id': ','.join(pubmed_ids),
        'retmode': 'xml'
    }
    print(f"[DEBUG] Fetching details for {len(pubmed_ids)} articles")
    response = requests.get(fetch_url, params=params)
    print(f"[DEBUG] efetch status code: {response.status_code}")
    return response.content

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
        abstract_elem = article.find('.//AbstractText')
        abstract = abstract_elem.text if (abstract_elem is not None and abstract_elem.text is not None) else ""
        
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
    results = []
    
    # Expanded list of protein-related keywords
    protein_keywords = ['protein', 'antigen', 'receptor', 'enzyme', 'kinase', 'factor', 'peptide', 
                       'ligand', 'antibody', 'cytokine', 'chemokine', 'growth factor', 'hormone']
    
    # Expanded list of cancer types
    cancer_types = [
        'breast cancer', 'lung cancer', 'prostate cancer', 'colorectal cancer', 'leukemia', 
        'melanoma', 'pancreatic cancer', 'ovarian cancer', 'carcinoma', 'sarcoma', 'glioblastoma',
        'bladder cancer', 'liver cancer', 'thyroid cancer', 'cervical cancer', 'endometrial cancer',
        'gastric cancer', 'esophageal cancer', 'head and neck cancer', 'brain cancer', 'lymphoma',
        'myeloma', 'renal cell carcinoma', 'testicular cancer', 'bone cancer'
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
        
        # Associate each protein with each cancer type found
        for protein in proteins_found:
            for cancer in cancers_found:
                results.append({
                    'protein': protein,
                    'cancer_type': cancer,
                    'pmid': article['pmid'],
                    'doi': article['doi']
                })
    
    print(f"[DEBUG] Found {len(results)} protein-cancer associations")
    return results

def main():
    # More specific query to get relevant articles
    query = "(cancer AND protein) AND (breast OR lung OR prostate OR colorectal OR leukemia OR pancreatic OR ovarian)"
    pubmed_ids = fetch_pubmed_articles(query, retmax=300)
    print(f"Found {len(pubmed_ids)} articles")
    
    # Add a delay to respect NCBI's rate limits
    time.sleep(1)
    
    xml_data = fetch_article_details(pubmed_ids)
    articles = parse_articles(xml_data)
    associations = extract_protein_cancer_associations(articles)
    
    # Deduplicate to unique protein-cancer pairs
    unique_associations = []
    seen = set()
    for assoc in associations:
        key = (assoc['protein'], assoc['cancer_type'])
        if key not in seen:
            seen.add(key)
            unique_associations.append(assoc)
    
    # Limit to 300 unique associations
    output_associations = unique_associations[:300]
    print(f"[DEBUG] Writing {len(output_associations)} unique associations to CSV")
    
    # Output to CSV
    with open('cancer_proteins.csv', 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['protein', 'cancer_type', 'pmid', 'doi']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for assoc in output_associations:
            writer.writerow(assoc)
    
    print("Data extraction complete. Output saved to cancer_proteins.csv")

if __name__ == "__main__":
    main()