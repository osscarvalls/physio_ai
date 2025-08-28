"""
Servicio directo de PubMed usando Entrez API
"""

import logging
import os
import json
from typing import List, Dict, Any, Optional
from Bio import Entrez
from datetime import datetime

from app.config.settings import settings
from qdrant_client import QdrantClient
from langchain_qdrant import QdrantVectorStore
from langchain_core.documents import Document
from langchain_openai import OpenAIEmbeddings

logger = logging.getLogger(__name__)


class PubMedService:
    """Servicio directo para consultar PubMed usando Entrez"""
    
    def __init__(self):
        """Inicializa el servicio de PubMed"""
        self._setup_entrez()
        self.docs_dir = "./docs"
    
    def _setup_entrez(self):
        """Configura Entrez para PubMed"""
        if settings.ENTREZ_EMAIL:
            Entrez.email = settings.ENTREZ_EMAIL
        else:
            logger.warning("ENTREZ_EMAIL no configurado. Usando email por defecto.")
            Entrez.email = "pysio_ai@example.com"
        
        Entrez.sleep_between_trials = 1
        Entrez.tool = "PySIO_AI/1.0"
    
    async def search_and_fetch(self, query: str, max_results: int = 5) -> List[Dict[str, Any]]:
        """
        Busca y obtiene artículos de PubMed
        
        Args:
            query: Consulta de búsqueda
            max_results: Número máximo de resultados
            
        Returns:
            Lista de artículos con detalles
        """
        try:
            logger.info(f"Buscando en PubMed: {query}")
            
            # Buscar PMIDs
            handle = Entrez.esearch(
                db='pubmed',
                sort='relevance',
                retmax=str(max_results),
                retmode='xml',
                term=query
            )
            results = Entrez.read(handle)
            handle.close()
            
            pmid_list = results.get('IdList', [])
            if not pmid_list:
                logger.warning(f"No se encontraron resultados para: {query}")
                return []
            
            # Obtener detalles de los artículos
            articles = []
            for pmid in pmid_list:
                try:
                    article = await self._fetch_article_details(pmid)
                    if article:
                        articles.append(article)
                        # Solo guardar en QDRANT si el artículo se obtuvo correctamente
                        self._save_article_to_qdrant(article)
                    else:
                        logger.warning(f"No se pudo obtener el artículo para PMID {pmid}")

                    # Pausa para no sobrecargar la API
                    import asyncio
                    await asyncio.sleep(0.1)
                    
                except Exception as e:
                    logger.error(f"Error procesando PMID {pmid}: {str(e)}")
                    continue
                
            
            logger.info(f"Procesados {len(articles)} artículos de {len(pmid_list)} PMIDs")
            return articles
            
        except Exception as e:
            logger.error(f"Error en búsqueda y obtención: {str(e)}")
            return []
    
    async def _fetch_article_details(self, pmid: str) -> Optional[Dict[str, Any]]:
        """Obtiene detalles de un artículo por PMID"""
        try:
            handle = Entrez.efetch(
                db="pubmed",
                id=pmid,
                rettype="xml",
                retmode="xml"
            )
            results = Entrez.read(handle)
            handle.close()
            
            if not results.get('PubmedArticle'):
                return None
            
            article = results['PubmedArticle'][0]
            return self._parse_article(article)
            
        except Exception as e:
            logger.error(f"Error obteniendo artículo {pmid}: {str(e)}")
            return None
    
    def _parse_article(self, article: Dict[str, Any]) -> Dict[str, Any]:
        """Parsea un artículo de PubMed"""
        try:
            pmid = str(article['MedlineCitation']['PMID'])
            title = article['MedlineCitation']['Article']['ArticleTitle']
            
            # Extraer abstract
            abstract = ""
            if 'Abstract' in article['MedlineCitation']['Article']:
                abstract_text = article['MedlineCitation']['Article']['Abstract'].get('AbstractText', [])
                if isinstance(abstract_text, list):
                    abstract = " ".join(abstract_text)
                else:
                    abstract = str(abstract_text)
            
            # Extraer autores
            authors = []
            if 'AuthorList' in article['MedlineCitation']['Article']:
                for author in article['MedlineCitation']['Article']['AuthorList']:
                    if 'LastName' in author and 'ForeName' in author:
                        authors.append(f"{author['ForeName']} {author['LastName']}")
            
            # Extraer fecha y journal
            pub_date = ""
            journal = ""
            if 'Journal' in article['MedlineCitation']['Article']:
                journal_info = article['MedlineCitation']['Article']['Journal']
                if 'Title' in journal_info:
                    journal = journal_info['Title']
                if 'PubDate' in journal_info:
                    pub_date_info = journal_info['PubDate']
                    if 'Year' in pub_date_info:
                        pub_date = pub_date_info['Year']
            
            # Extraer palabras clave
            keywords = []
            if 'MeshHeadingList' in article['MedlineCitation']:
                for mesh in article['MedlineCitation']['MeshHeadingList']:
                    if 'DescriptorName' in mesh:
                        keywords.append(mesh['DescriptorName'])
            
            # Extraer full text
            full_text = ""
            if 'FullJournalName' in article['MedlineCitation']['Article']['Journal']:
                full_text = article['MedlineCitation']['Article']['Journal']['FullJournalName']
            
            return {
                'pmid': pmid,
                'title': title,
                'abstract': abstract,
                'authors': authors,
                'publication_date': pub_date,
                'journal': journal,
                'keywords': keywords,
                'full_text': full_text
            }
            
        except Exception as e:
            logger.error(f"Error parseando artículo: {str(e)}")
            return {}
    
    def _save_article_to_qdrant(self, article: Dict[str, Any]):
        """Guarda un artículo en Qdrant"""
        try:
            # Validar que el artículo tenga los campos mínimos necesarios
            if not article or not isinstance(article, dict):
                logger.error("Artículo inválido o None")
                return
                
            if 'pmid' not in article:
                logger.error("Artículo sin PMID")
                return
                
            # Asegurar que page_content no esté vacío
            page_content = article.get("abstract", "")
            if not page_content:
                page_content = article.get("title", "Sin contenido disponible")
            
            # Crear cliente de Qdrant
            client = QdrantClient(
                url=settings.QDRANT_URL,
                api_key=settings.QDRANT_API_KEY
            )
            # Crear vectorstore de Qdrant
            collection_name = settings.QDRANT_COLLECTION_NAME
            vectorstore = QdrantVectorStore(
                client=client,
                collection_name=collection_name,
                embedding=OpenAIEmbeddings(
                    api_key=settings.OPENAI_API_KEY
                )
            )
            # Crear documento simple sin ID
            doc = Document(
                page_content=page_content,
                metadata={
                    "title": article.get("title", ""),
                    "pmid": article['pmid'],
                    "keywords": article.get("keywords", ""),
                    "authors": article.get("authors", ""),
                    "publication_date": article.get("publication_date", ""),
                    "journal": article.get("journal", "")
                }
            )
            
            # Debug: verificar que el documento se creó correctamente
            logger.info(f"Documento creado: PMID={article['pmid']}, contenido={len(page_content)} chars")
            logger.info(f"Metadatos: {doc.metadata}")

            # Intentar agregar usando langchain primero
            try:
                vectorstore.add_documents([doc])
                logger.info(f"Artículo {article['pmid']} guardado exitosamente en QDRANT usando langchain")
            except Exception as add_error:
                logger.warning(f"Langchain falló, intentando con cliente directo: {str(add_error)}")
                
                # Fallback: usar cliente directo de QDRANT
                try:
                    # Crear embeddings
                    embeddings = OpenAIEmbeddings(api_key=settings.OPENAI_API_KEY)
                    vector = embeddings.embed_query(page_content)
                    
                    # Insertar directamente en QDRANT
                    client.upsert(
                        collection_name=collection_name,
                        points=[{
                            "id": article['pmid'],
                            "vector": vector,
                            "payload": {
                                "page_content": page_content,
                                "title": article.get("title", ""),
                                "pmid": article['pmid'],
                                "keywords": article.get("keywords", ""),
                                "authors": article.get("authors", ""),
                                "publication_date": article.get("publication_date", ""),
                                "journal": article.get("journal", "")
                            }
                        }]
                    )
                    logger.info(f"Artículo {article['pmid']} guardado exitosamente en QDRANT usando cliente directo")
                    
                except Exception as direct_error:
                    logger.error(f"Error con cliente directo: {str(direct_error)}")

        except Exception as e:
            logger.error(f"Error guardando el artículo {article['pmid']} en QDRANT: {str(e)}")
    