import os
import json

from api.engine.utils import set_openai, get_prompts

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
        self.pubmed_search_template = ChatPromptTemplate.from_template(self.prompts["pubmed_search"])

    
    def generate_pubmed_search_entries(self,symptoms):
        ## Use LLM to generate search entries to query pubmed API
        print(self.prompts["pubmed_search"])
        pm_search_llm = self.pubmed_search_template | self.llm
        search_entries = pm_search_llm.invoke({'symptoms':symptoms})
        search_entries = search_entries.content.replace('[','').replace(']','').replace("'","").split(', ')
        print("Search entries: ", search_entries)
        return search_entries

    def search_pubmed(self, entry):
        ## Search pubmed API
        Entrez.email = 'oscarvallslozano@gmail.com'
        handle = Entrez.esearch(db='pubmed',
                                sort='relevance',
                                retmax='10',
                                retmode='xml',
                                term=entry)
        results = Entrez.read(handle)
        return results
    
    def fetch_papers(self,id_list):
        ## Fetch papers from pubmed API
        ids = ','.join(id_list)
        Entrez.email = 'oscarvallslozano@gmail.com'
        handle = Entrez.efetch(db='pubmed',
                            retmode='xml',
                            id=ids)
        results = Entrez.read(handle)
        return results
    
    def store_docs(self,papers):
        ## Store all papers fetched from pubmed API in local directory
        for paper in papers['PubmedArticle']:
            try:
                paper_id = paper['MedlineCitation']['PMID']
                paper_title = paper['MedlineCitation']['Article']['ArticleTitle']
                abstract_text = paper['MedlineCitation']['Article']['Abstract']['AbstractText']
                if not os.path.isfile(f'./docs/{paper_id}.txt'):
                    with open(file=f'./docs/{paper_id}.txt',mode='w') as f:
                        f.write(f"Title: {paper_title}\n")
                        f.write(f"Abstract: {abstract_text}")
                else:
                    print(f"Document {paper_id} already exists")
            except Exception as e:
                print(f"Error storing doc {paper_id}: {e}")
                continue
    
    def retrieve_evidence(self,search_entries):
        ## Generate search entries and search pubmed API for each search entry and store the results in local docs
        for entry in search_entries:
            results = self.search_pubmed(entry)
            print(entry)
            print(results)
            id_list = results['IdList']
            papers = self.fetch_papers(id_list)
            self.store_docs(papers)

    def __call__(self,symptoms):
        ## Callable function to execute evidence retrieval and context generation
        ## STEPS:
        ## 1. Generate search entries
        search_entries = self.generate_pubmed_search_entries(symptoms)
        ## 2. Search pubmed API and store results in local docs
        self.retrieve_evidence(search_entries)
        ## 3. Load and split docs
        loader = DirectoryLoader("./docs",glob="**/*.txt")
        docs = loader.load()
        text_splitter = RecursiveCharacterTextSplitter(chunk_size=1000, chunk_overlap=200, add_start_index=True) 
        all_splits = text_splitter.split_documents(docs)
        ## 4. Create vector store
        vectorstore = Chroma.from_documents(documents=all_splits, embedding=OpenAIEmbeddings())
        ## 5. Create retriever
        retriever = vectorstore.as_retriever(search_type="similarity", search_kwargs={"k": 10})
        ## 6. Return retriever
        return retriever