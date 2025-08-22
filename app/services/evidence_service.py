"""
Servicio para recuperar evidencia médica de PubMed
"""

import logging
import os
from typing import List, Dict, Any, Optional
from Bio import Entrez
from langchain.schema import Document
from langchain_community.document_loaders import DirectoryLoader
from langchain_text_splitters import RecursiveCharacterTextSplitter
from langchain_chroma import Chroma
from langchain_openai import OpenAIEmbeddings

from app.config.settings import settings
from app.models.diagnosis import MedicalEvidence

logger = logging.getLogger(__name__)


class EvidenceService:
    """Servicio para recuperar y gestionar evidencia médica"""
    
    def __init__(self):
        """Inicializa el servicio de evidencia"""
        self._setup_entrez()
        self._setup_vector_store()
        self.text_splitter = RecursiveCharacterTextSplitter(
            chunk_size=1000,
            chunk_overlap=200
        )
    
    def _setup_entrez(self):
        """Configura Entrez para PubMed"""
        if settings.ENTREZ_EMAIL:
            Entrez.email = settings.ENTREZ_EMAIL
        else:
            logger.warning("ENTREZ_EMAIL no configurado. Usando email por defecto.")
            Entrez.email = "pysio_ai@example.com"
    
    def _setup_vector_store(self):
        """Configura la base de datos vectorial"""
        try:
            self.embeddings = OpenAIEmbeddings(
                api_key=settings.OPENAI_API_KEY
            )
            self.vector_store = Chroma(
                persist_directory=settings.CHROMA_PERSIST_DIR,
                embedding_function=self.embeddings
            )
        except Exception as e:
            logger.error(f"Error configurando vector store: {str(e)}")
            self.vector_store = None
    
    async def search_pubmed(self, query: str, max_results: int = 10) -> List[Dict[str, Any]]:
        """
        Busca artículos en PubMed
        
        Args:
            query: Consulta de búsqueda
            max_results: Número máximo de resultados
            
        Returns:
            Lista de resultados de PubMed
        """
        try:
            logger.info(f"Buscando en PubMed: {query}")
            
            # Realizar búsqueda
            handle = Entrez.esearch(
                db='pubmed',
                sort='relevance',
                retmax=str(max_results),
                retmode='xml',
                term=query
            )
            results = Entrez.read(handle)
            handle.close()
            
            return results.get('IdList', [])
            
        except Exception as e:
            logger.error(f"Error buscando en PubMed: {str(e)}")
            return []
    
    async def fetch_papers(self, pmid_list: List[str]) -> List[Dict[str, Any]]:
        """
        Obtiene detalles de artículos por PMID
        
        Args:
            pmid_list: Lista de IDs de PubMed
            
        Returns:
            Lista de artículos con detalles
        """
        if not pmid_list:
            return []
        
        try:
            logger.info(f"Obteniendo {len(pmid_list)} artículos de PubMed")
            
            handle = Entrez.efetch(
                db="pubmed",
                id=pmid_list,
                rettype="xml",
                retmode="xml"
            )
            results = Entrez.read(handle)
            handle.close()
            
            return results.get('PubmedArticle', [])
            
        except Exception as e:
            logger.error(f"Error obteniendo artículos: {str(e)}")
            return []
    
    async def process_papers(self, papers: List[Dict[str, Any]]) -> List[MedicalEvidence]:
        """
        Procesa artículos de PubMed y los convierte a modelos de evidencia
        
        Args:
            papers: Lista de artículos de PubMed
            
        Returns:
            Lista de evidencia médica procesada
        """
        evidence_list = []
        
        for paper in papers:
            try:
                # Extraer información del artículo
                pmid = paper['MedlineCitation']['PMID']
                title = paper['MedlineCitation']['Article']['ArticleTitle']
                
                # Manejar abstract (puede no existir)
                abstract = ""
                if 'Abstract' in paper['MedlineCitation']['Article']:
                    abstract_text = paper['MedlineCitation']['Article']['Abstract'].get('AbstractText', [])
                    if isinstance(abstract_text, list):
                        abstract = " ".join(abstract_text)
                    else:
                        abstract = str(abstract_text)
                
                # Crear modelo de evidencia
                evidence = MedicalEvidence(
                    pmid=pmid,
                    title=title,
                    abstract=abstract,
                    relevance_score=0.0  # Se puede implementar scoring real
                )
                
                evidence_list.append(evidence)
                
            except Exception as e:
                logger.error(f"Error procesando artículo {pmid}: {str(e)}")
                continue
        
        return evidence_list
    
    async def store_papers_locally(self, papers: List[Dict[str, Any]]) -> None:
        """
        Almacena artículos localmente para uso posterior
        
        Args:
            papers: Lista de artículos de PubMed
        """
        if not os.path.exists(settings.DOCS_DIR):
            os.makedirs(settings.DOCS_DIR)
        
        for paper in papers:
            try:
                pmid = paper['MedlineCitation']['PMID']
                title = paper['MedlineCitation']['Article']['ArticleTitle']
                
                # Obtener abstract
                abstract = ""
                if 'Abstract' in paper['MedlineCitation']['Article']:
                    abstract_text = paper['MedlineCitation']['Article']['Abstract'].get('AbstractText', [])
                    if isinstance(abstract_text, list):
                        abstract = " ".join(abstract_text)
                    else:
                        abstract = str(abstract_text)
                
                # Guardar archivo
                file_path = os.path.join(settings.DOCS_DIR, f"{pmid}.txt")
                if not os.path.exists(file_path):
                    with open(file_path, 'w', encoding='utf-8') as f:
                        f.write(f"Title: {title}\n")
                        f.write(f"Abstract: {abstract}\n")
                        f.write(f"PMID: {pmid}\n")
                    
                    logger.info(f"Artículo {pmid} guardado localmente")
                else:
                    logger.debug(f"Artículo {pmid} ya existe localmente")
                    
            except Exception as e:
                logger.error(f"Error guardando artículo {pmid}: {str(e)}")
                continue
    
    async def retrieve_relevant_evidence(
        self, 
        symptoms: str, 
        max_results: int = 5
    ) -> List[Document]:
        """
        Recupera evidencia médica relevante para los síntomas
        
        Args:
            symptoms: Síntomas del paciente
            max_results: Número máximo de documentos a retornar
            
        Returns:
            Lista de documentos relevantes
        """
        try:
            if not self.vector_store:
                logger.warning("Vector store no disponible, usando búsqueda local")
                return await self._search_local_docs(symptoms)
            
            # Búsqueda semántica en vector store
            results = self.vector_store.similarity_search(
                symptoms, 
                k=max_results
            )
            
            return results
            
        except Exception as e:
            logger.error(f"Error recuperando evidencia: {str(e)}")
            return await self._search_local_docs(symptoms)
    
    async def _search_local_docs(self, symptoms: str) -> List[Document]:
        """Búsqueda fallback en documentos locales"""
        try:
            if not os.path.exists(settings.DOCS_DIR):
                return []
            
            # Cargar documentos locales
            loader = DirectoryLoader(
                settings.DOCS_DIR,
                glob="**/*.txt",
                show_progress=True
            )
            documents = loader.load()
            
            # Dividir en chunks
            split_docs = self.text_splitter.split_documents(documents)
            
            # Retornar primeros documentos (simulación de relevancia)
            return split_docs[:5]
            
        except Exception as e:
            logger.error(f"Error en búsqueda local: {str(e)}")
            return []
    
    async def update_vector_store(self) -> None:
        """Actualiza la base de datos vectorial con documentos locales"""
        try:
            if not self.vector_store:
                logger.warning("Vector store no disponible")
                return
            
            # Cargar documentos locales
            if not os.path.exists(settings.DOCS_DIR):
                logger.info("Directorio de documentos no existe")
                return
            
            loader = DirectoryLoader(
                settings.DOCS_DIR,
                glob="**/*.txt",
                show_progress=True
            )
            documents = loader.load()
            
            if not documents:
                logger.info("No hay documentos para procesar")
                return
            
            # Dividir documentos
            split_docs = self.text_splitter.split_documents(documents)
            
            # Actualizar vector store
            self.vector_store.add_documents(split_docs)
            self.vector_store.persist()
            
            logger.info(f"Vector store actualizado con {len(split_docs)} documentos")
            
        except Exception as e:
            logger.error(f"Error actualizando vector store: {str(e)}")
