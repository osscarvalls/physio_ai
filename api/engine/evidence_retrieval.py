import os
import json

from utils import set_openai, get_prompts

from langchain_openai import ChatOpenAI
from langchain.prompts import ChatPromptTemplate
from langchain_community.document_loaders import DirectoryLoader
from langchain_text_splitters import RecursiveCharacterTextSplitter
from langchain_chroma import Chroma
from langchain_openai import OpenAIEmbeddings

from Bio import Entrez

class EvidenceRetrievalEngine():
    def __init__(self) -> None:
        ## Load environment variables
        set_openai()
        self.prompts = get_prompts()
        self.llm = ChatOpenAI(model='gpt-4o-mini',temperature=0)       
        self.pubmed_search_template = ChatPromptTemplate.from_template(self.prompts["pubmed_search_prompt"])

    
    def generate_pubmed_search_entries(self,symptoms):
        pm_search_llm = self.pubmed_search_template | self.llm
        search_entries = pm_search_llm.invoke({'symptoms':symptoms})
        return search_entries.content

    def search_pubmed(self, entry):
        Entrez.email = 'oscarvallslozano@gmail.com'
        handle = Entrez.esearch(db='pubmed',
                                sort='relevance',
                                retmax='10',
                                retmode='xml',
                                term=entry)
        results = Entrez.read(handle)
        return results
    
    def fetch_papers(self,id_list):
        ids = ','.join(id_list)
        Entrez.email = 'oscarvallslozano@gmail.com'
        handle = Entrez.efetch(db='pubmed',
                            retmode='xml',
                            id=ids)
        results = Entrez.read(handle)
        return results
    
    def store_docs(self,papers):
        for paper in papers['PubmedArticle']:
            paper_id = paper['MedlineCitation']['PMID']
            paper_title = paper['MedlineCitation']['Article']['ArticleTitle']
            abstract_text = paper['MedlineCitation']['Article']['Abstract']['AbstractText']
            if not os.path.isfile(f'./docs/{paper_id}.txt'):
                with open(file=f'./docs/{paper_id}.txt',mode='w') as f:
                    f.write(f"Title: {paper_title}\n")
                    f.write(f"Abstract: {abstract_text}")
            else:
                continue

    def retrieve_evidence(self,search_entries):
        for entry in search_entries:
            results = self.search_pubmed(entry)
            id_list = results['IdList']
            papers = self.fetch_papers(id_list)
            self.store_docs(papers)

    def generate_vector_store(self):
        loader = DirectoryLoader("./docs",glob="**/*.txt")
        docs = loader.load()
        text_splitter = RecursiveCharacterTextSplitter(chunk_size=1000, chunk_overlap=200, add_start_index=True) 
        all_splits = text_splitter.split_documents(docs)
        vectorstore = Chroma.from_documents(documents=all_splits, embedding=OpenAIEmbeddings())