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
                        await self._save_article_locally(article)
                    
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
            pmid = article['MedlineCitation']['PMID']
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
            
            return {
                'pmid': pmid,
                'title': title,
                'abstract': abstract,
                'authors': authors,
                'publication_date': pub_date,
                'journal': journal,
                'keywords': keywords
            }
            
        except Exception as e:
            logger.error(f"Error parseando artículo: {str(e)}")
            return {}
    
    async def _save_article_locally(self, article: Dict[str, Any]) -> bool:
        """Guarda un artículo localmente"""
        try:
            os.makedirs(self.docs_dir, exist_ok=True)
            
            pmid = article['pmid']
            safe_title = "".join(c for c in article['title'] if c.isalnum() or c in (' ', '-', '_')).rstrip()
            safe_title = safe_title[:50]
            
            filename = f"{pmid}_{safe_title}.txt"
            filepath = os.path.join(self.docs_dir, filename)
            
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(f"PMID: {pmid}\n")
                f.write(f"Title: {article['title']}\n")
                f.write(f"Journal: {article['journal']}\n")
                f.write(f"Publication Date: {article['publication_date']}\n")
                f.write(f"Authors: {', '.join(article['authors'])}\n")
                f.write(f"Keywords: {', '.join(article['keywords'])}\n")
                f.write("\n" + "="*80 + "\n\n")
                f.write("ABSTRACT:\n")
                f.write(article['abstract'])
            
            logger.info(f"Artículo {pmid} guardado en {filepath}")
            return True
            
        except Exception as e:
            logger.error(f"Error guardando artículo {article.get('pmid', 'unknown')}: {str(e)}")
            return False
    
    async def get_local_articles(self) -> List[Dict[str, Any]]:
        """Obtiene artículos guardados localmente"""
        try:
            if not os.path.exists(self.docs_dir):
                return []
            
            articles = []
            for filename in os.listdir(self.docs_dir):
                if filename.endswith('.txt'):
                    filepath = os.path.join(self.docs_dir, filename)
                    
                    with open(filepath, 'r', encoding='utf-8') as f:
                        content = f.read()
                    
                    # Extraer PMID del contenido
                    pmid = ""
                    for line in content.split('\n'):
                        if line.startswith('PMID:'):
                            pmid = line.split(':', 1)[1].strip()
                            break
                    
                    articles.append({
                        'pmid': pmid,
                        'content': content,
                        'filepath': filepath
                    })
            
            return articles
            
        except Exception as e:
            logger.error(f"Error obteniendo artículos locales: {str(e)}")
            return []
