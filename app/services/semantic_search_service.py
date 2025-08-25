"""
Servicio de búsqueda semántica para evidencia médica
"""

import logging
import os
from typing import List, Dict, Any
from langchain_openai import OpenAIEmbeddings
from langchain_chroma import Chroma
from langchain_text_splitters import RecursiveCharacterTextSplitter

from app.config.settings import settings
from app.services.pubmed_service import PubMedService

logger = logging.getLogger(__name__)


class SemanticSearchService:
    """Servicio de búsqueda semántica para evidencia médica"""
    
    def __init__(self):
        """Inicializa el servicio de búsqueda semántica"""
        self.pubmed_service = PubMedService()
        self._setup_components()
    
    def _setup_components(self):
        """Configura los componentes de LangChain"""
        try:
            # Configurar embeddings
            self.embeddings = OpenAIEmbeddings(
                api_key=settings.OPENAI_API_KEY
            )
            
            # Configurar vector store
            self.vector_store = Chroma(
                persist_directory=settings.CHROMA_PERSIST_DIR,
                embedding_function=self.embeddings
            )
            
            # Configurar text splitter
            self.text_splitter = RecursiveCharacterTextSplitter(
                chunk_size=1000,
                chunk_overlap=200
            )
            
            logger.info("Componentes de búsqueda semántica configurados")
            
        except Exception as e:
            logger.error(f"Error configurando componentes: {str(e)}")
            self.embeddings = None
            self.vector_store = None
            self.text_splitter = None
    
    async def ingest_pubmed_articles(self, query: str, max_results: int = 5) -> bool:
        """
        Ingresa artículos de PubMed en la base de datos vectorial
        
        Args:
            query: Consulta de búsqueda para PubMed
            max_results: Número máximo de artículos a obtener
            
        Returns:
            True si se ingresaron correctamente
        """
        try:
            logger.info(f"Ingestando artículos de PubMed para: {query}")
            
            # Obtener artículos de PubMed
            articles = await self.pubmed_service.search_and_fetch(query, max_results)
            
            if not articles:
                logger.warning("No se encontraron artículos para ingestar")
                return False
            
            # Preparar documentos para LangChain
            documents = []
            for article in articles:
                # Crear contenido del documento
                content = f"Title: {article['title']}\n"
                content += f"Abstract: {article['abstract']}\n"
                content += f"Journal: {article['journal']}\n"
                content += f"Keywords: {', '.join(article['keywords'])}"
                
                # Crear metadatos
                metadata = {
                    "pmid": article['pmid'],
                    "title": article['title'],
                    "journal": article['journal'],
                    "publication_date": article['publication_date'],
                    "authors": ", ".join(article['authors']),
                    "source": "pubmed",
                    "query": query
                }
                
                # Crear documento de LangChain
                from langchain_core.documents import Document
                doc = Document(
                    page_content=content,
                    metadata=metadata
                )
                documents.append(doc)
            
            # Dividir documentos en chunks
            if self.text_splitter:
                split_docs = self.text_splitter.split_documents(documents)
                logger.info(f"Documentos divididos en {len(split_docs)} chunks")
            else:
                split_docs = documents
                logger.warning("Text splitter no disponible, usando documentos completos")
            
            # Añadir a la base de datos vectorial
            if self.vector_store:
                self.vector_store.add_documents(split_docs)
                self.vector_store.persist()
                logger.info(f"Se ingresaron {len(split_docs)} chunks en la base de datos vectorial")
                return True
            else:
                logger.error("Vector store no disponible")
                return False
                
        except Exception as e:
            logger.error(f"Error ingestando artículos: {str(e)}")
            return False
    
    async def search_medical_evidence(self, symptoms: str, max_results: int = 5) -> List[Dict[str, Any]]:
        """
        Busca evidencia médica relevante para los síntomas
        
        Args:
            symptoms: Síntomas del paciente
            max_results: Número máximo de resultados
            
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
            
            # Formatear resultados
            formatted_results = []
            for i, doc in enumerate(results):
                formatted_results.append({
                    'content': doc.page_content,
                    'metadata': doc.metadata,
                    'rank': i + 1,
                    'source': 'vector_store'
                })
            
            logger.info(f"Búsqueda semántica completada, {len(formatted_results)} resultados")
            return formatted_results
            
        except Exception as e:
            logger.error(f"Error en búsqueda semántica: {str(e)}")
            return await self._search_local_docs(symptoms)
    
    async def _search_local_docs(self, symptoms: str) -> List[Dict[str, Any]]:
        """Búsqueda fallback en documentos locales"""
        try:
            # Obtener artículos locales
            local_articles = await self.pubmed_service.get_local_articles()
            
            if not local_articles:
                return []
            
            # Búsqueda simple por palabras clave
            relevant_articles = []
            symptoms_lower = symptoms.lower()
            
            for article in local_articles:
                content_lower = article['content'].lower()
                
                # Calcular relevancia simple
                relevance_score = 0
                symptom_words = symptoms_lower.split()
                
                for word in symptom_words:
                    if len(word) > 3:  # Ignorar palabras muy cortas
                        if word in content_lower:
                            relevance_score += 1
                
                if relevance_score > 0:
                    relevant_articles.append({
                        'content': article['content'],
                        'metadata': {
                            'pmid': article['pmid'],
                            'source': 'local',
                            'relevance_score': relevance_score
                        },
                        'rank': len(relevant_articles) + 1,
                        'source': 'local_search'
                    })
            
            # Ordenar por relevancia
            relevant_articles.sort(key=lambda x: x['metadata']['relevance_score'], reverse=True)
            
            logger.info(f"Búsqueda local completada, {len(relevant_articles)} resultados")
            return relevant_articles[:5]  # Limitar a 5 resultados
            
        except Exception as e:
            logger.error(f"Error en búsqueda local: {str(e)}")
            return []
    
    async def get_vector_store_stats(self) -> Dict[str, Any]:
        """
        Obtiene estadísticas de la base de datos vectorial
        
        Returns:
            Dict con estadísticas
        """
        try:
            if not self.vector_store:
                return {"error": "Vector store no disponible"}
            
            # Obtener información de la colección
            collection = self.vector_store._collection
            count = collection.count()
            
            return {
                "total_documents": count,
                "collection_name": collection.name,
                "status": "active"
            }
            
        except Exception as e:
            logger.error(f"Error obteniendo estadísticas: {str(e)}")
            return {"error": str(e)}
