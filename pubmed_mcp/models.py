"""
Modelos Pydantic para artículos de PubMed
"""

from typing import List, Optional
from pydantic import BaseModel, Field


class PubMedArticle(BaseModel):
    """Modelo para un artículo de PubMed"""
    
    pmid: str = Field(..., description="PMID (PubMed ID) único del artículo")
    title: str = Field(..., description="Título del artículo")
    abstract: str = Field(default="", description="Resumen/abstract del artículo")
    authors: List[str] = Field(default_factory=list, description="Lista de autores del artículo")
    publication_date: str = Field(default="", description="Fecha de publicación (año)")
    journal: str = Field(default="", description="Nombre de la revista/journal")
    keywords: List[str] = Field(default_factory=list, description="Palabras clave MeSH del artículo")
    full_text: str = Field(default="", description="Nombre completo de la revista")
    
    class Config:
        """Configuración del modelo"""
        json_schema_extra = {
            "example": {
                "pmid": "33190432",
                "title": "Physiotherapy and rehabilitation applications in lipedema management: A literature review.",
                "abstract": "Lipedema is a chronic and progressive disease of adipose tissue...",
                "authors": ["M Esmer", "F J Schingale", "D Unal"],
                "publication_date": "2019",
                "journal": "Lymphology",
                "keywords": ["Combined Modality Therapy", "Exercise Therapy", "Humans"],
                "full_text": "Lymphology"
            }
        }


class PubMedSearchResponse(BaseModel):
    """Modelo para la respuesta de búsqueda en PubMed"""
    
    articles: List[PubMedArticle] = Field(..., description="Lista de artículos encontrados")
    total_results: int = Field(..., description="Número total de resultados")
    query: str = Field(..., description="Consulta de búsqueda realizada")
    
    class Config:
        """Configuración del modelo"""
        json_schema_extra = {
            "example": {
                "articles": [
                    {
                        "pmid": "33190432",
                        "title": "Physiotherapy and rehabilitation applications...",
                        "abstract": "Lipedema is a chronic...",
                        "authors": ["M Esmer", "F J Schingale"],
                        "publication_date": "2019",
                        "journal": "Lymphology",
                        "keywords": ["Combined Modality Therapy"],
                        "full_text": "Lymphology"
                    }
                ],
                "total_results": 1,
                "query": "physiotherapy rehabilitation"
            }
        }
