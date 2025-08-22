"""
Configuración centralizada de la aplicación PySIO AI
"""

import os
from typing import Optional
from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    """Configuración de la aplicación usando Pydantic Settings"""
    
    # Configuración del servidor
    HOST: str = "0.0.0.0"
    PORT: int = 8000
    DEBUG: bool = False
    
    # Configuración de OpenAI
    OPENAI_API_KEY: str
    OPENAI_ORGANIZATION_ID: Optional[str] = None
    OPENAI_MODEL: str = "gpt-4o-mini"
    OPENAI_TEMPERATURE: float = 0.0
    
    # Configuración de PubMed
    ENTREZ_EMAIL: Optional[str] = None
    
    # Configuración de archivos y directorios
    STATIC_DIR: str = "templates/static"
    TEMPLATES_DIR: str = "templates"
    DOCS_DIR: str = "docs"
    
    # Configuración de la base de datos de documentos
    CHROMA_PERSIST_DIR: str = "./chroma_db"
    
    class Config:
        env_file = ".env"
        case_sensitive = True


# Instancia global de configuración
settings = Settings()
