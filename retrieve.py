from Bio import Entrez

# Set your email (NCBI requires this)
Entrez.email = "your_email@example.com"

# Define OR-based search term for "bicep curls", "exercise angles", "time under tension", etc.
query = "bicep curls OR exercise angles OR time under tension OR muscle hypertrophy"
handle = Entrez.esearch(db="pubmed", term=query, retmax=100)  # Adjust 'retmax' for more results
record = Entrez.read(handle)
handle.close()

# Get PubMed IDs (PMIDs) of matching articles
pmids = record["IdList"]
print("PubMed IDs:", pmids)

# Fetch article details (abstracts)
def fetch_article_details(pmids):
    """Retrieve article details from PubMed given a list of PMIDs."""
    handle = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract", retmode="text")
    articles = handle.read()
    handle.close()
    return articles

# Fetch and display articles
article_data = fetch_article_details(pmids)
print(article_data)
